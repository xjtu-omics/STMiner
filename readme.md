# Usage

### import package
```python
import scanpy as sc
import squidpy as sq

from util import *
from Algorithm.graph import *
from Algorithm.distribution import *
```

### load data
```python
h5ad_path = 'F://Rep11_MOB_ST.h5ad'
adata = sc.read_h5ad(h5ad_path)
```
### preprocess
```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
```
### fit GMM 
```python
# get high variable gene 
gene_list = adata.var[adata.var['highly_variable']==True].index
# fit gmms
gmm_dict = fit_gmms(adata, gene_list, n_comp=10, thread=4)
```
### build distance matrix
```python
arr = build_distance_array(gmm_dict)
```
### clustering
```python
cluster()
```
