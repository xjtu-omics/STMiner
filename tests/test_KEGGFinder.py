import pandas as pd
import pytest

from STMiner import KEGGFinder, KEGGFinderError


class FakeServices:
    TIMEOUT = 30
    url = "http://rest.kegg.jp"


class FakeKEGG:
    def __init__(self):
        self.services = FakeServices()
        self.responses = {}

    def get(self, entry_id):
        return self.responses.get(entry_id, 404)

    @staticmethod
    def parse(raw):
        if raw == "generic":
            return {"ENTRY": "map00010 Pathway", "COMPOUND": {"C00031": "Glucose"}}
        if raw == "genes":
            return {
                "ENTRY": "hsa00010 Pathway",
                "GENE": {
                    "7157": "TP53; tumor protein p53",
                    "1": "A1BG",
                },
            }
        return {}

    @staticmethod
    def get_pathway_by_gene(gene, organism):
        if gene == "7157" and organism == "hsa":
            return {"hsa04115": "p53 signaling pathway"}
        return None

    @staticmethod
    def parse_kgml_pathway(pathway):
        return {"entries": [pathway], "relations": []}

    @staticmethod
    def pathway2sif(pathway, uniprot=True):
        return [(pathway, "interacts_with", str(uniprot))]


class FakeQuickGO:
    def __init__(self):
        self.services = FakeServices()
        self.calls = []

    def Annotation(self, **kwargs):
        self.calls.append(kwargs)
        page = kwargs["page"]
        return {
            "results": [
                {
                    "geneProductId": "UniProtKB:P04637",
                    "symbol": "TP53",
                    "goId": f"GO:{page:07d}",
                    "goName": f"term {page}",
                    "goAspect": "biological_process",
                }
            ],
            "pageInfo": {"current": page, "total": 2},
        }


class FakeUniProt:
    def __init__(self):
        self.services = FakeServices()

    @staticmethod
    def search(query, **kwargs):
        return "Entry\tGene Names\nP04637\tTP53"

    @staticmethod
    def mapping(**kwargs):
        return {
            "results": [{"from": "P04637", "to": "hsa:7157"}],
            "failedIds": ["BAD"],
        }


@pytest.fixture
def finder():
    kegg = FakeKEGG()
    kegg.responses.update({"map00010": "generic", "hsa00010": "genes"})
    return KEGGFinder(
        timeout=12,
        kegg=kegg,
        quickgo=FakeQuickGO(),
        uniprot=FakeUniProt(),
    )


def test_find_validates_input_and_clears_stale_result(finder):
    with pytest.raises(ValueError, match="non-empty"):
        finder.find("  ")

    finder.find("hsa00010")
    with pytest.raises(LookupError, match="was not found"):
        finder.find("hsa99999")
    assert finder.result is None
    assert finder.entry_id is None


def test_generic_entry_and_gene_table_are_stable(finder):
    result = finder.find("map00010")
    assert result["ENTRY"] == "map00010 Pathway"
    assert finder.get_gene_dataframe().empty
    assert finder.get_section_dataframe("compound").to_dict("records") == [
        {"id": "C00031", "name": "Glucose", "info": ""}
    ]
    with pytest.raises(LookupError, match="no GENE section"):
        finder.get_gene_dataframe(strict=True)


def test_gene_table_handles_missing_semicolon(finder):
    finder.find("hsa00010")
    frame = finder.get_gene_dataframe()
    assert frame.to_dict("records") == [
        {"id": "7157", "symbol": "TP53", "info": "tumor protein p53"},
        {"id": "1", "symbol": "A1BG", "info": ""},
    ]


def test_common_kegg_interfaces(finder):
    pathways = finder.get_pathways_by_gene("hsa:7157", "hsa")
    assert pathways.to_dict("records") == [
        {"pathway_id": "hsa04115", "name": "p53 signaling pathway"}
    ]
    assert finder.get_pathway_network("hsa04115")["entries"] == ["hsa04115"]
    assert finder.pathway_to_sif("hsa04115", uniprot=False) == [
        ("hsa04115", "interacts_with", "False")
    ]
    with pytest.raises(ValueError, match="organism prefix"):
        finder.get_pathways_by_gene("mmu:7157", "hsa")


def test_quickgo_annotations_are_paginated_and_normalised(finder):
    frame = finder.get_go_annotations(["P04637", "P04637"], taxon_id=9606, aspect="BP")
    assert frame.shape[0] == 2
    assert frame["goId"].tolist() == ["GO:0000001", "GO:0000002"]
    assert len(finder.quickgo.calls) == 2
    assert finder.quickgo.calls[0]["aspect"] == "P"
    assert finder.quickgo.services.TIMEOUT == 12


def test_uniprot_search_and_mapping(finder):
    search = finder.search_uniprot("gene:TP53")
    assert search.to_dict("records") == [{"Entry": "P04637", "Gene Names": "TP53"}]

    mapping = finder.map_identifiers(["P04637", "BAD"])
    assert mapping.to_dict("records") == [{"from": "P04637", "to": "hsa:7157"}]
    assert mapping.attrs["failed_ids"] == ["BAD"]


def test_enrichment_delegates_to_scanpy_gprofiler(monkeypatch, finder):
    import scanpy

    result_frame = pd.DataFrame(
        {
            "source": ["KEGG"],
            "native": ["KEGG:04115"],
            "name": ["p53 signaling pathway"],
            "p_value": [0.001],
            "intersections": [["TP53"]],
        }
    )
    captured = {}

    def fake_enrich(query, *, org, gprofiler_kwargs):
        captured.update(
            query=query,
            org=org,
            gprofiler_kwargs=gprofiler_kwargs,
        )
        return result_frame

    monkeypatch.setattr(scanpy.queries, "enrich", fake_enrich)
    result = finder.enrich_kegg(["TP53", "TP53", "BRCA1"])

    assert result.equals(result_frame)
    assert captured["query"] == ["TP53", "BRCA1"]
    assert captured["org"] == "hsapiens"
    assert captured["gprofiler_kwargs"]["sources"] == ["KEGG"]
    assert captured["gprofiler_kwargs"]["no_evidences"] is False
    assert finder.enrichment_result.equals(result_frame)


def test_enrichment_rejects_invalid_arguments(finder):
    with pytest.raises(ValueError, match="at least one identifier"):
        finder.enrich_gene_set([])
    with pytest.raises(ValueError, match="correction_method"):
        finder.enrich_gene_set(["TP53"], correction_method="BH")
    with pytest.raises(TypeError, match="all_results"):
        finder.enrich_gene_set(["TP53"], all_results=1)


def test_unexpected_remote_responses_raise_stable_error(finder):
    finder.quickgo.Annotation = lambda **kwargs: "bad response"
    with pytest.raises(KEGGFinderError, match="unexpected response"):
        finder.get_go_annotations(["P04637"])


@pytest.fixture
def enrichment_result():
    return pd.DataFrame(
        {
            "source": ["GO:BP", "GO:BP", "KEGG", "KEGG"],
            "native": ["GO:1", "GO:2", "KEGG:1", "KEGG:2"],
            "name": [
                "DNA repair",
                "cell cycle checkpoint signaling",
                "p53 signaling pathway",
                "Homologous recombination",
            ],
            "p_value": [0.001, 0.01, 0.0001, 0.02],
            "intersection_size": [4, 3, 5, 2],
            "query_size": [10, 10, 10, 10],
            "query": ["damage", "damage", "damage", "stemness"],
        }
    )


def test_plot_enrichment_creates_nature_style_panels_and_export(
    finder, enrichment_result, tmp_path
):
    original = enrichment_result.copy(deep=True)
    output = tmp_path / "enrichment.svg"

    fig, axes = finder.plot_enrichment(
        enrichment_result,
        top_n=2,
        save_path=output,
    )

    assert output.exists() and output.stat().st_size > 0
    assert len(axes) == 2
    assert [axis.get_title(loc="left") for axis in axes] == [
        "GO biological process",
        "KEGG pathway",
    ]
    assert fig.get_size_inches()[0] == pytest.approx(7.2)
    pd.testing.assert_frame_equal(enrichment_result, original)
    from matplotlib import pyplot as plt

    plt.close(fig)


def test_plot_enrichment_uses_latest_result_and_validates_data(
    finder, enrichment_result
):
    import matplotlib.pyplot as plt

    finder.enrichment_result = enrichment_result
    fig, axes = finder.plot_enrichment(sources="KEGG", top_n=1)
    assert len(axes) == 1
    plt.close(fig)

    finder.enrichment_result = None
    with pytest.raises(RuntimeError, match="Run an enrichment method"):
        finder.plot_enrichment()
    with pytest.raises(ValueError, match="required columns"):
        finder.plot_enrichment(pd.DataFrame({"name": ["term"]}))


def test_plot_enrichment_uses_gene_ratio_for_one_gene_set(finder, enrichment_result):
    from matplotlib import pyplot as plt

    single_query = enrichment_result.drop(columns="query")
    fig, axes = finder.plot_enrichment(single_query, sources="KEGG", top_n=2)

    assert len(axes) == 1
    assert axes[0].get_xlabel() == "Gene ratio"
    plt.close(fig)
