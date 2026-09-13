from __future__ import annotations

from collections.abc import Iterable, Mapping
from io import StringIO
from typing import Any

import pandas as pd
from bioservices import KEGG, QuickGO, UniProt


class KEGGFinderError(RuntimeError):
    """Raised when a remote annotation service cannot complete a request."""


class KEGGFinder:
    """Query KEGG, QuickGO, UniProt, and g:Profiler through stable wrappers."""

    _GO_ASPECTS = {
        "P": "P",
        "BP": "P",
        "BIOLOGICAL_PROCESS": "P",
        "F": "F",
        "MF": "F",
        "MOLECULAR_FUNCTION": "F",
        "C": "C",
        "CC": "C",
        "CELLULAR_COMPONENT": "C",
    }
    _GO_COLUMNS = [
        "geneProductId",
        "symbol",
        "goId",
        "goName",
        "goAspect",
        "qualifier",
        "goEvidence",
        "evidenceCode",
        "reference",
        "taxonId",
        "assignedBy",
        "date",
    ]

    def __init__(
        self,
        timeout: int = 30,
        verbose: bool = False,
        cache: bool = False,
        *,
        kegg: Any | None = None,
        quickgo: Any | None = None,
        uniprot: Any | None = None,
    ):
        if not isinstance(timeout, int) or isinstance(timeout, bool) or timeout <= 0:
            raise ValueError("timeout must be a positive integer.")

        self.timeout = timeout
        self.verbose = verbose
        self.cache = cache
        self._kegg = kegg
        self._quickgo = quickgo
        self._uniprot = uniprot
        self.result: dict[str, Any] | None = None
        self.entry_id: str | None = None
        self.enrichment_result: pd.DataFrame | None = None

        for client in (kegg, quickgo, uniprot):
            self._set_timeout(client)

    @property
    def kegg(self):
        if self._kegg is None:
            self._kegg = KEGG(verbose=self.verbose, cache=self.cache)
            services = getattr(self._kegg, "services", None)
            if services is not None and hasattr(services, "url"):
                services.url = "https://rest.kegg.jp"
            self._set_timeout(self._kegg)
        return self._kegg

    @property
    def quickgo(self):
        if self._quickgo is None:
            self._quickgo = QuickGO(verbose=self.verbose, cache=self.cache)
            self._set_timeout(self._quickgo)
        return self._quickgo

    @property
    def uniprot(self):
        if self._uniprot is None:
            self._uniprot = UniProt(verbose=self.verbose, cache=self.cache)
            self._set_timeout(self._uniprot)
        return self._uniprot

    def _set_timeout(self, client: Any | None) -> None:
        services = getattr(client, "services", None)
        if services is not None and hasattr(services, "TIMEOUT"):
            services.TIMEOUT = self.timeout

    @staticmethod
    def _validate_text(value: Any, name: str) -> str:
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{name} must be a non-empty string.")
        return value.strip()

    @staticmethod
    def _normalise_values(values: Iterable[Any], name: str) -> list[str]:
        if isinstance(values, (str, bytes)):
            raise TypeError(f"{name} must be an iterable of identifiers, not a string.")
        try:
            normalised = [str(value).strip() for value in values]
        except TypeError as exc:
            raise TypeError(f"{name} must be an iterable of identifiers.") from exc
        normalised = list(dict.fromkeys(value for value in normalised if value))
        if not normalised:
            raise ValueError(f"{name} must contain at least one identifier.")
        return normalised

    @staticmethod
    def _raise_service_error(service: str, operation: str, exc: Exception) -> None:
        raise KEGGFinderError(
            f"{service} request failed while {operation}: {exc}"
        ) from exc

    def find(self, pathway: str) -> dict[str, Any]:
        """Retrieve and parse one KEGG entry while preserving the legacy API."""
        return self.get_entry(pathway)

    def get_entry(self, entry_id: str) -> dict[str, Any]:
        """Retrieve any text-based KEGG entry, including reference entries."""
        entry_id = self._validate_text(entry_id, "entry_id")
        self.result = None
        self.entry_id = None

        try:
            raw = self.kegg.get(entry_id)
        except Exception as exc:
            self._raise_service_error("KEGG", f"retrieving {entry_id!r}", exc)

        if not isinstance(raw, str) or not raw.strip():
            status = raw if isinstance(raw, int) else "empty response"
            raise LookupError(f"KEGG entry {entry_id!r} was not found ({status}).")

        try:
            parsed = self.kegg.parse(raw)
        except Exception as exc:
            self._raise_service_error("KEGG", f"parsing {entry_id!r}", exc)

        if not isinstance(parsed, dict) or not parsed:
            raise KEGGFinderError(
                f"KEGG returned an unparseable entry for {entry_id!r}."
            )

        self.result = parsed
        self.entry_id = entry_id
        return parsed

    def get_section_dataframe(self, section: str) -> pd.DataFrame:
        """Convert a mapping-like KEGG entry section into a stable table."""
        section = self._validate_text(section, "section").upper()
        if self.result is None:
            raise RuntimeError("Call find() or get_entry() before reading a section.")

        values = self.result.get(section)
        columns = ["id", "name", "info"]
        if values is None:
            return pd.DataFrame(columns=columns)
        if not isinstance(values, Mapping):
            raise TypeError(
                f"KEGG section {section!r} is not mapping-like and cannot be tabulated."
            )

        records = []
        for identifier, description in values.items():
            name, separator, info = str(description).partition(";")
            records.append(
                {
                    "id": str(identifier),
                    "name": name.strip(),
                    "info": info.strip() if separator else "",
                }
            )
        return pd.DataFrame.from_records(records, columns=columns)

    def get_gene_dataframe(self, strict: bool = False) -> pd.DataFrame:
        """Return the GENE section using the historical id/symbol/info columns."""
        if not isinstance(strict, bool):
            raise TypeError("strict must be a boolean.")
        frame = self.get_section_dataframe("GENE")
        if frame.empty and strict:
            entry = self.entry_id or "current KEGG entry"
            raise LookupError(f"KEGG entry {entry!r} has no GENE section.")
        return frame.rename(columns={"name": "symbol"})

    def get_pathways_by_gene(self, gene: str, organism: str) -> pd.DataFrame:
        """Return KEGG pathways containing one organism-specific gene."""
        gene = self._validate_text(gene, "gene")
        organism = self._validate_text(organism, "organism")
        if ":" in gene:
            gene_organism, gene = gene.split(":", maxsplit=1)
            if gene_organism != organism:
                raise ValueError(
                    "The organism prefix in gene must match the organism argument."
                )

        try:
            pathways = self.kegg.get_pathway_by_gene(gene, organism)
        except Exception as exc:
            self._raise_service_error(
                "KEGG", f"finding pathways for {organism}:{gene}", exc
            )

        if pathways is None:
            return pd.DataFrame(columns=["pathway_id", "name"])
        if not isinstance(pathways, Mapping):
            raise KEGGFinderError("KEGG returned an unexpected pathway result.")
        records = [
            {"pathway_id": str(pathway_id), "name": str(name)}
            for pathway_id, name in pathways.items()
        ]
        return pd.DataFrame.from_records(records, columns=["pathway_id", "name"])

    def get_pathway_network(self, pathway: str) -> dict[str, Any]:
        """Return parsed KGML entries and relations for a KEGG pathway."""
        pathway = self._validate_text(pathway, "pathway")
        try:
            network = self.kegg.parse_kgml_pathway(pathway)
        except Exception as exc:
            self._raise_service_error("KEGG", f"parsing KGML for {pathway!r}", exc)
        if not isinstance(network, dict):
            raise KEGGFinderError("KEGG returned an unexpected KGML result.")
        return network

    def pathway_to_sif(self, pathway: str, uniprot: bool = False) -> list[Any]:
        """Convert a KEGG pathway to the bioservices SIF representation."""
        pathway = self._validate_text(pathway, "pathway")
        if not isinstance(uniprot, bool):
            raise TypeError("uniprot must be a boolean.")
        try:
            result = self.kegg.pathway2sif(pathway, uniprot=uniprot)
        except Exception as exc:
            self._raise_service_error("KEGG", f"converting {pathway!r} to SIF", exc)
        if result is None:
            return []
        if not isinstance(result, list):
            raise KEGGFinderError("KEGG returned an unexpected SIF result.")
        return result

    def get_go_annotations(
        self,
        gene_products: Iterable[Any],
        *,
        taxon_id: int | str | None = None,
        aspect: str | None = None,
        page_size: int = 100,
        max_pages: int | None = None,
        batch_size: int = 50,
    ) -> pd.DataFrame:
        """Retrieve paginated QuickGO annotations for gene product identifiers."""
        gene_products = self._normalise_values(gene_products, "gene_products")
        if taxon_id is not None:
            if isinstance(taxon_id, bool) or not isinstance(taxon_id, (int, str)):
                raise TypeError("taxon_id must be an integer, string, or None.")
            taxon_id = str(taxon_id).strip()
            if not taxon_id.isdigit() or int(taxon_id) <= 0:
                raise ValueError("taxon_id must be a positive taxonomy identifier.")
        if (
            not isinstance(page_size, int)
            or isinstance(page_size, bool)
            or not 1 <= page_size <= 100
        ):
            raise ValueError("page_size must be an integer between 1 and 100.")
        if (
            not isinstance(batch_size, int)
            or isinstance(batch_size, bool)
            or batch_size <= 0
        ):
            raise ValueError("batch_size must be a positive integer.")
        if max_pages is not None and (
            not isinstance(max_pages, int)
            or isinstance(max_pages, bool)
            or max_pages <= 0
        ):
            raise ValueError("max_pages must be a positive integer or None.")

        quickgo_aspect = None
        if aspect is not None:
            key = self._validate_text(aspect, "aspect").upper()
            if key not in self._GO_ASPECTS:
                raise ValueError("aspect must be one of BP, MF, CC, P, F, or C.")
            quickgo_aspect = self._GO_ASPECTS[key]

        records: list[dict[str, Any]] = []
        for start in range(0, len(gene_products), batch_size):
            batch = gene_products[start : start + batch_size]
            page = 1
            while True:
                parameters = {
                    "geneProductId": ",".join(batch),
                    "includeFields": "goName,taxonName",
                    "limit": page_size,
                    "page": page,
                }
                if taxon_id is not None:
                    parameters["taxonId"] = str(taxon_id)
                if quickgo_aspect is not None:
                    parameters["aspect"] = quickgo_aspect

                try:
                    response = self.quickgo.Annotation(**parameters)
                except Exception as exc:
                    self._raise_service_error(
                        "QuickGO", "retrieving gene annotations", exc
                    )
                if not isinstance(response, dict):
                    raise KEGGFinderError("QuickGO returned an unexpected response.")

                page_records = response.get("results", [])
                if not isinstance(page_records, list):
                    raise KEGGFinderError("QuickGO returned malformed annotation rows.")
                records.extend(page_records)

                page_info = response.get("pageInfo") or {}
                total_pages = page_info.get("total", 1)
                if not isinstance(total_pages, int):
                    total_pages = 1
                if page >= total_pages or (max_pages is not None and page >= max_pages):
                    break
                page += 1

        if not records:
            return pd.DataFrame(columns=self._GO_COLUMNS)
        frame = pd.DataFrame.from_records(records).reindex(columns=self._GO_COLUMNS)
        return frame.drop_duplicates(ignore_index=True)

    def search_uniprot(
        self,
        query: str,
        *,
        columns: str = "accession,id,gene_names,organism_name",
        limit: int | None = None,
    ) -> pd.DataFrame:
        """Run a UniProt search and parse its TSV response."""
        query = self._validate_text(query, "query")
        columns = self._validate_text(columns, "columns")
        if limit is not None and (
            not isinstance(limit, int) or isinstance(limit, bool) or limit <= 0
        ):
            raise ValueError("limit must be a positive integer or None.")
        try:
            response = self.uniprot.search(
                query,
                frmt="tsv",
                columns=columns,
                limit=limit,
                progress=False,
            )
        except Exception as exc:
            self._raise_service_error("UniProt", f"searching for {query!r}", exc)
        if response is None or response == "":
            return pd.DataFrame()
        if not isinstance(response, str):
            raise KEGGFinderError("UniProt returned an unexpected search response.")
        try:
            return pd.read_csv(StringIO(response), sep="\t")
        except Exception as exc:
            raise KEGGFinderError("UniProt returned malformed TSV data.") from exc

    def map_identifiers(
        self,
        identifiers: Iterable[Any],
        *,
        source: str = "UniProtKB_AC-ID",
        target: str = "KEGG",
    ) -> pd.DataFrame:
        """Map identifiers with the UniProt ID mapping service."""
        identifiers = self._normalise_values(identifiers, "identifiers")
        source = self._validate_text(source, "source")
        target = self._validate_text(target, "target")
        try:
            response = self.uniprot.mapping(
                fr=source,
                to=target,
                query=",".join(identifiers),
                progress=False,
            )
        except Exception as exc:
            self._raise_service_error(
                "UniProt", f"mapping {source!r} to {target!r}", exc
            )
        if not isinstance(response, dict):
            raise KEGGFinderError("UniProt returned an unexpected mapping response.")
        results = response.get("results", [])
        if not isinstance(results, list):
            raise KEGGFinderError("UniProt returned malformed mapping rows.")
        frame = (
            pd.json_normalize(results)
            if results
            else pd.DataFrame(columns=["from", "to"])
        )
        frame.attrs["failed_ids"] = response.get("failedIds", [])
        return frame

    def enrich_gene_set(
        self,
        genes: Iterable[Any] | Mapping[str, Iterable[Any]],
        *,
        organism: str = "hsapiens",
        sources: Iterable[str] = ("GO:BP", "GO:MF", "GO:CC", "KEGG"),
        background: Iterable[Any] | None = None,
        user_threshold: float = 0.05,
        correction_method: str = "g_SCS",
        all_results: bool = False,
        ordered: bool = False,
        no_iea: bool = False,
    ) -> pd.DataFrame:
        """Run GO and KEGG enrichment through Scanpy's g:Profiler wrapper."""
        if isinstance(genes, Mapping):
            if not genes:
                raise ValueError("genes must contain at least one named gene set.")
            query: list[str] | dict[str, list[str]] = {
                self._validate_text(name, "gene set name"): self._normalise_values(
                    values, f"genes[{name!r}]"
                )
                for name, values in genes.items()
            }
        else:
            query = self._normalise_values(genes, "genes")

        organism = self._validate_text(organism, "organism")
        sources = self._normalise_values(sources, "sources")
        if not isinstance(user_threshold, (int, float)) or isinstance(
            user_threshold, bool
        ):
            raise TypeError("user_threshold must be numeric.")
        if not 0 < user_threshold <= 1:
            raise ValueError("user_threshold must be greater than 0 and at most 1.")
        if correction_method not in {"g_SCS", "bonferroni", "fdr"}:
            raise ValueError(
                "correction_method must be 'g_SCS', 'bonferroni', or 'fdr'."
            )
        for value, name in (
            (all_results, "all_results"),
            (ordered, "ordered"),
            (no_iea, "no_iea"),
        ):
            if not isinstance(value, bool):
                raise TypeError(f"{name} must be a boolean.")

        kwargs: dict[str, Any] = {
            "sources": sources,
            "user_threshold": user_threshold,
            "significance_threshold_method": correction_method,
            "all_results": all_results,
            "ordered": ordered,
            "no_iea": no_iea,
            "no_evidences": False,
        }
        if background is not None:
            kwargs["background"] = self._normalise_values(background, "background")
            kwargs["domain_scope"] = "custom"

        try:
            from scanpy import queries as scanpy_queries

            result = scanpy_queries.enrich(
                query,
                org=organism,
                gprofiler_kwargs=kwargs,
            )
        except ImportError as exc:
            raise KEGGFinderError(
                "Gene-set enrichment requires gprofiler-official. "
                "Install the project dependencies before calling this method."
            ) from exc
        except Exception as exc:
            self._raise_service_error("g:Profiler", "running gene-set enrichment", exc)

        if not isinstance(result, pd.DataFrame):
            raise KEGGFinderError("g:Profiler returned an unexpected result.")
        self.enrichment_result = result.copy()
        return result

    def enrich_go(
        self,
        genes: Iterable[Any] | Mapping[str, Iterable[Any]],
        **kwargs: Any,
    ) -> pd.DataFrame:
        """Run GO BP, MF, and CC enrichment through g:Profiler."""
        kwargs["sources"] = ("GO:BP", "GO:MF", "GO:CC")
        return self.enrich_gene_set(genes, **kwargs)

    def enrich_kegg(
        self,
        genes: Iterable[Any] | Mapping[str, Iterable[Any]],
        **kwargs: Any,
    ) -> pd.DataFrame:
        """Run KEGG enrichment through g:Profiler."""
        kwargs["sources"] = ("KEGG",)
        return self.enrich_gene_set(genes, **kwargs)

    def cluster_gene_set(
        self,
        genes: Iterable[Any] | Mapping[str, Iterable[Any]],
        **kwargs: Any,
    ) -> pd.DataFrame:
        """Return functional GO/KEGG clusters for one or more gene sets."""
        return self.enrich_gene_set(genes, **kwargs)

    def plot_enrichment(
        self,
        result: pd.DataFrame | None = None,
        **kwargs: Any,
    ):
        """Plot a gene-set enrichment result with the Nature-style module."""
        if result is None:
            if self.enrichment_result is None:
                raise RuntimeError(
                    "Run an enrichment method or provide result before plotting."
                )
            result = self.enrichment_result

        from STMiner.Plot.enrichment import plot_enrichment

        return plot_enrichment(result, **kwargs)
