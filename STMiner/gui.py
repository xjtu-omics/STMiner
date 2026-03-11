import importlib
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from STMiner.SPFinder import SPFinder


def build_spfinder_gui(spfinder: "SPFinder"):
    if importlib.util.find_spec("gradio") is None:
        raise ImportError(
            "Gradio is required for the GUI. Please install it with `pip install gradio`."
        )

    import gradio as gr

    def _check_adata_loaded():
        if spfinder.adata is None:
            raise ValueError("Please load spatial transcriptomics data first.")

    def load_data(file, file_path, amplification, bin_size, merge_bin):
        try:
            path = None
            if file is not None:
                path = file.name
            elif file_path:
                path = file_path
            if path is None:
                return "Please provide an h5ad file path or upload a file."
            amplification = int(amplification)
            bin_size = int(bin_size)
            spfinder.read_h5ad(
                file=path,
                amplification=amplification,
                bin_size=bin_size,
                merge_bin=merge_bin,
            )
            obs_num, var_num = spfinder.adata.n_obs, spfinder.adata.n_vars
            return f"Loaded data from {path}. Spots: {obs_num}, Genes: {var_num}."
        except Exception as exc:
            return f"Failed to load data: {exc}"

    def find_spatial_hvg(min_cells, min_genes, vmax, thread):
        try:
            _check_adata_loaded()
            min_cells = int(min_cells)
            min_genes = int(min_genes)
            vmax = int(vmax)
            thread = int(thread)
            spfinder.get_genes_csr_array(
                min_cells=min_cells,
                min_genes=min_genes,
                vmax=vmax,
            )
            spfinder.spatial_high_variable_genes(vmax=vmax, thread=thread)
            preview = spfinder.global_distance.head(20)
            return "Computed spatially variable genes.", preview
        except Exception as exc:
            return f"Failed to compute: {exc}", None

    def fit_patterns(n_top_genes, n_comp, min_cells, gene_list_text):
        try:
            _check_adata_loaded()
            n_top_genes = int(n_top_genes)
            n_comp = int(n_comp)
            min_cells = int(min_cells)
            genes = None
            if gene_list_text:
                genes = [g.strip() for g in gene_list_text.split("\n") if g.strip()]
            spfinder.fit_pattern(
                n_top_genes=n_top_genes,
                n_comp=n_comp,
                min_cells=min_cells,
                gene_list=genes,
            )
            return "Finished fitting GMM patterns."
        except Exception as exc:
            return f"Failed to fit patterns: {exc}"

    def cluster_genes(method, n_clusters, mds_components, use_highly_variable_gene, n_top_genes, gene_list_text):
        try:
            _check_adata_loaded()
            n_clusters = int(n_clusters)
            mds_components = int(mds_components)
            n_top_genes = int(n_top_genes)
            genes = None
            if gene_list_text:
                genes = [g.strip() for g in gene_list_text.split("\n") if g.strip()]
            spfinder.build_distance_array(method=method, gene_list=genes)
            spfinder.cluster_gene(
                n_clusters=n_clusters,
                mds_components=mds_components,
                use_highly_variable_gene=use_highly_variable_gene,
                n_top_genes=n_top_genes,
            )
            preview = None
            if spfinder.genes_labels is not None:
                preview = spfinder.genes_labels.head(20)
            return "Clustering completed.", preview
        except Exception as exc:
            return f"Failed to cluster genes: {exc}", None

    with gr.Blocks(title="STMiner SPFinder GUI") as demo:
        gr.Markdown(
            """
            # STMiner.SPFinder GUI

            load spatial data, find spatially variable genes,
            fit Gaussian Mixture Models, build distance arrays, and cluster genes.
            """
        )

        with gr.Tab("Step.1-Load data"):
            gr.Markdown(
                "Upload or provide an h5ad path, then configure binning options to load the dataset."
            )
            with gr.Row():
                file_input = gr.File(label="Upload h5ad", file_types=[".h5ad"])
                path_input = gr.Textbox(label="Or enter h5ad path", placeholder="/path/to/file.h5ad")
            with gr.Row():
                amplification = gr.Number(label="Amplification", value=1, precision=0)
                bin_size = gr.Number(label="Bin size", value=1, precision=0)
                merge_bin = gr.Checkbox(label="Merge bin", value=False)
            load_btn = gr.Button("Load dataset")
            load_status = gr.Textbox(label="Status")
            load_btn.click(
                load_data,
                inputs=[file_input, path_input, amplification, bin_size, merge_bin],
                outputs=load_status,
            )

        with gr.Tab("Step.2-Spatial HVG"):
            gr.Markdown("Compute spatially variable genes following the README workflow.")
            with gr.Row():
                min_cells = gr.Number(label="Min cells", value=100, precision=0)
                min_genes = gr.Number(label="Min genes", value=10, precision=0)
                vmax = gr.Number(label="Vmax percentile", value=100, precision=0)
                thread = gr.Number(label="Threads", value=1, precision=0)
            hvg_btn = gr.Button("Find spatial HVGs")
            hvg_status = gr.Textbox(label="Status")
            hvg_preview = gr.DataFrame(label="Top spatially variable genes", interactive=False)
            hvg_btn.click(
                find_spatial_hvg,
                inputs=[min_cells, min_genes, vmax, thread],
                outputs=[hvg_status, hvg_preview],
            )

        with gr.Tab("Step.3-Fit patterns"):
            gr.Markdown("Fit GMM patterns for genes.")
            with gr.Row():
                n_top_genes = gr.Number(label="Top highly variable genes (-1 for all)", value=2000, precision=0)
                n_comp = gr.Number(label="GMM components", value=20, precision=0)
                min_cells_fit = gr.Number(label="Min cells", value=20, precision=0)
            gene_list_fit = gr.Textbox(
                label="Optional gene list (one per line)",
                placeholder="GeneA\nGeneB",
            )
            fit_btn = gr.Button("Fit patterns")
            fit_status = gr.Textbox(label="Status")
            fit_btn.click(
                fit_patterns,
                inputs=[n_top_genes, n_comp, min_cells_fit, gene_list_fit],
                outputs=fit_status,
            )

        with gr.Tab("Step.4-Distance & clustering"):
            gr.Markdown("Build distance arrays and cluster genes into spatial patterns.")
            with gr.Row():
                method = gr.Radio(
                    label="Distance method",
                    choices=["gmm", "mse", "cs", "ot"],
                    value="gmm",
                )
                n_clusters = gr.Number(label="Number of clusters", value=6, precision=0)
                mds_components = gr.Number(label="MDS components", value=20, precision=0)
            use_hvg = gr.Checkbox(label="Use highly variable genes for clustering", value=False)
            n_top_cluster = gr.Number(label="Top genes when using HVG", value=500, precision=0)
            gene_list_cluster = gr.Textbox(
                label="Optional gene list for distance array (one per line)",
                placeholder="GeneA\nGeneB",
            )
            cluster_btn = gr.Button("Build & cluster")
            cluster_status = gr.Textbox(label="Status")
            cluster_preview = gr.DataFrame(label="Gene labels", interactive=False)
            cluster_btn.click(
                cluster_genes,
                inputs=[
                    method,
                    n_clusters,
                    mds_components,
                    use_hvg,
                    n_top_cluster,
                    gene_list_cluster,
                ],
                outputs=[cluster_status, cluster_preview],
            )

    return demo


def launch_spfinder_gui(spfinder: "SPFinder"):
    demo = build_spfinder_gui(spfinder)
    return demo.launch()
