# Usage

## Quick start by example

### import package

```python
from SPFinder import SPFinder
```

### Load data

```python
spf = SPFinder()
spf.read_10x(file='F://Rep11_MOB_ST.h5ad', amplification=1000, bin_size=80)
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
spf.cluster(n_clusters=6)
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
plot_heatmap(label=0, vmax=95)
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
- read_10x
- read_gem
- merge_bin
#### preprocess
- fit_pattern
- normalize
- log1p
#### fit model
- fit_pattern
#### build distance array & clustering
- cluster
#### visualization
