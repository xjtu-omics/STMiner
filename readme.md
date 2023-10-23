<div align=center><img src="./pic/logo.png" height = "200"/></div>

# Introduction
Spatial transcriptomics revolutionizes transcriptomics by incorporating positional information. However, an emergency problem is to find out the gene expression pattern which can reveal the special region in tissue and find out the genes only expression in those regions.
![STMiner](./pic/fig1.png)

Here we propose “STMiner” based on the Gaussian mixture model to solve this problem. STMiner is a bottom-up methodology algorithm. It is initiated by fitting a parametric model of gene spatial distributions and constructing a distance array between them utilizing the Hellinger distance. Genes are clustered, thereby recognizing spatial co-expression patterns across distinct gene classes.
Please visit [STMiner Documents](https://stminerdoc.readthedocs.io/en/latest/Introduction/Introduction.html) for
details.

## Quick start by example

### import package

```python
from STMiner.SPFinder import SPFinder
```

### Load data

```python
spf = SPFinder()
spf.read_h5ad(file='F://Rep11_MOB_ST.h5ad', amplification=1000, bin_size=80)
```

### Preprocess

```python
spf.normalize()
```

Optional:

```python
spf.log1p()
```

### Fit GMM

```python
spf.fit_pattern(n_top_genes=100, n_comp=10)
```

Each GMM model has 10 components.

### build distance matrix & clustering

```python
spf.cluster_gene(n_clusters=6)
```

### Result & Visualization

The result are stored in **genes_labels**:

```python
spf.genes_labels
```

The output looks like the following:

|    | gene_id        | labels |
|----|----------------|--------|
| 0  | Cldn5          | 2      |
| 1  | Fyco1          | 2      |
| 2  | Pmepa1         | 2      |
| 3  | Arhgap5        | 0      |
| 4  | Apc            | 5      |
| .. | ...            | ...    |
| 95 | Cyp2a5         | 0      |
| 96 | X5730403I07Rik | 0      |
| 97 | Ltbp2          | 2      |
| 98 | Rbp4           | 4      |
| 99 | Hist1h1e       | 4      |

To visualize the patterns by heatmap:

```python
spf.plot_pattern(vmax=95)
```

To visualize the genes expression heatmap by labels:

```python
plot_genes(label=0, vmax=95)
```

## API

### Attribute of SPFinder

| Attribute           | Type         | Description                        |
|---------------------|--------------|------------------------------------|
| adata               | Anndata      | Anndata for loaded spatial data    |
| genes_patterns      | dict         | GMM model for each gene            |
| genes_distance_aray | pd.DataFrame | Distance between each GMM          |
| genes_labels        | pd.DataFrame | Gene name and their pattern labels |

### Methods of SPFinder

#### load data

- read_h5ad
- read_gem
- merge_bin

#### preprocess

- fit_pattern
- normalize
- log1p

#### fit model

- fit_pattern

#### build distance array & clustering

- cluster_gene

#### visualization
