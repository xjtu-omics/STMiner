# Usage

### import package

```python
from SPFinder import SPFinder
```

### load data

```python
spf = SPFinder()
spf.read_10x(file='F://Rep11_MOB_ST.h5ad', amplification=1000, bin_size=80)
```

### preprocess

```python
spf.normalize()
```

Optional:

```python
spf.log1p()
```

### fit GMM

```python
spf.fit_pattern(100, n_comp=10)
```

### build distance matrix & clustering

```python
spf.cluster(6)
```

### Visualization

```python
spf.plot_pattern()
```

| Attribute           | Type         | Description                        |
|---------------------|--------------|------------------------------------|
| adata               | Anndata      | Anndata for loaded spatial data    |
| genes_patterns      | dict         | GMM model for each gene            |
| genes_distance_aray | pd.DataFrame | distance between each GMM          |
| genes_labels        | pd.DataFrame | gene name and their pattern labels |
