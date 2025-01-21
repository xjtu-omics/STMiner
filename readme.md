![Static Badge](https://img.shields.io/badge/License-GPL3.0-blue)
![PyPI - Version](https://img.shields.io/pypi/v/STMiner)
![GitHub repo size](https://img.shields.io/github/repo-size/xjtu-omics/STMiner)
![PyPI - Downloads](https://img.shields.io/pypi/dm/STMiner)
![Static Badge](https://img.shields.io/badge/3.10-green?logo=python&label=Python&labelColor=yellow)
![Static Badge](https://img.shields.io/badge/Linux-blue?logo=Linux&logoColor=white)
![Static Badge](https://img.shields.io/badge/Windows-blue?logo=Windows&logoColor=white)
![Static Badge](https://img.shields.io/badge/macos-blue?logo=apple&logoColor=white)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14678661.svg)](https://doi.org/10.5281/zenodo.14678661)
[![Documentation Status](https://readthedocs.org/projects/stminerdoc/badge/?version=latest)](https://stminerdoc.readthedocs.io/en/latest/?badge=latest) 
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/xjtu-omics/STMiner)
![GitHub Repo stars](https://img.shields.io/github/stars/xjtu-omics/STMiner)

<div align=center><img src="./pic/logo.png" height = "200"/></div>

# 👩‍🏫 Introduction
## Why STMiner?

ST data presents challenges such as uneven cell density distribution, low sampling rates, and complex spatial structures. Traditional spot-based analysis strategies struggle to effectively address these issues. STMiner explores ST data by leveraging the spatial distribution of genes, thus avoiding the biases that these conditions can introduce into the results.<br>

Most importantly, STMiner offers seamless integration with [Anndata](https://github.com/scverse/anndata)/[Scanpy](https://github.com/scverse/scanpy) and can be **easily installed** via [PyPI](https://pypi.org/project/STMiner/).<br>

<div align=center><img src="./pic/why.png" height = "500"/></div>

## Method detail
Here we propose “**STMiner**”. The three key steps of analyzing ST data in STMiner are depicted. 

<div align=center><img src="./pic/fig1.png" height = "400"/></div>

(**Left top**) STMiner first utilizes Gaussian Mixture Models (GMMs) to represent the spatial distribution of each gene and the overall spatial distribution. (**Left bottom**) STMiner then identifies spatially variable genes by calculating the cost that transfers the overall spatial distribution to gene spatial distribution. Genes with high costs exhibit significant spatial variation, meaning their expression patterns differ considerably across different regions of the tissue. The distance array is built between SVGs in the same way, genes with similar spatial structures have a low cost to transport to each other, and vice versa. (**Right**) The distance array is embedded into a low-dimensional space by Multidimensional Scaling, allowing for clustering genes with similar spatial expression patterns into distinct functional gene sets and getting their spatial structure. 



# 🚀 Quick start by example
**Please visit [STMiner Documents](https://stminerdoc.readthedocs.io/en/latest/Introduction/Introduction.html) for installation and detail usage.**

## import package

```python
from STMiner import SPFinder
```

## Load data

You can download the demo dataset from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4838133), or you can also download them from [STMOMICS](https://db.cngb.org/stomics/datasets/STDS0000086).
STMiner can read spatial transcriptome data in various formats, such as **gem**, **bmk**, and **h5ad** (see [STMiner Documents](https://stminerdoc.readthedocs.io/en/latest/Introduction/Introduction.html)).   
We recommend using the **h5ad** format, as it is currently the most widely used and supported by most algorithms and software in the spatial transcriptomics field.

```python
sp = SPFinder()
file_path = 'Path/to/your/h5ad/file'
sp.read_h5ad(file=file_path, bin_size=1)
```
## Find spatial high variable genes

```python
sp.get_genes_csr_array(min_cells=500, log1p=False)
sp.spatial_high_variable_genes()
```
 - The parameter **min_cells** was used to filter genes that are too sparse to generate a reliable spatial distribution.
 - The parameter **log1p** was used to avoid extreme values affecting the results. For most open-source h5ad files, log1p has already been executed, so the default value here is False.
 - You can perform STMiner in your interested gene sets. Use parameter **gene_list** to input the gene list to STMiner. Then, STMiner will only calculate the given gene set of the dataset.

You can check the distance of each gene by:

```python
sp.global_distance
```

| Gene  | Distance | z-score |
| ----- | -------- | ------- |
| geneA | 9998     | 5.5     |
| geneB | 9994     | 5.4     |
| ...   | ...      | 5.3     |
| geneC | 8724     | 5.2     |

The first column is the gene name, and the second column is the difference between the spatial distribution of the gene and the background.</br>
A larger difference indicates a more pronounced spatial pattern of the gene.

## Preprocess and Fit GMM

```python
sp.fit_pattern(n_comp=20, gene_list=list(sp.global_distance[:1000]['Gene']))
```

**n_comp=20** means each GMM model has 20 components.

## Build distance matrix & clustering

```python
# This step calculates the distance between genes' spatial distributions.
sp.build_distance_array()
# Dimensionality reduction and clustering.
sp.cluster_gene(n_clusters=6, mds_components=20) 
```

## Result & Visualization

The result is stored in **genes_labels**:

```python
sp.genes_labels
```

The output looks like the following:

|     | gene_id        | labels |
| --- | -------------- | ------ |
| 0   | Cldn5          | 2      |
| 1   | Fyco1          | 2      |
| 2   | Pmepa1         | 2      |
| 3   | Arhgap5        | 0      |
| 4   | Apc            | 5      |
| ..  | ...            | ...    |
| 95  | Cyp2a5         | 0      |
| 96  | X5730403I07Rik | 0      |
| 97  | Ltbp2          | 2      |
| 98  | Rbp4           | 4      |
| 99  | Hist1h1e       | 4      |

### Visualize the distance array:

```python
import seaborn as sns
sns.clustermap(sp.genes_distance_array)
```
<div align=center><img src="./pic/heatmap.png" width = "400"/></div>

### To visualize the patterns:
**Note**: A image path for ***image_path*** is needed if you want to show background image. In this example, you can download the processed image [here](https://github.com/xjtu-omics/STMiner/blob/main/pic/demo_img.png). Anyway, ***image_path*** is **optional**, not providing background images has no impact on the calculation results.

```python
sp.get_pattern_array(vote_rate=0.3)
img_path = 'path/to/downloaded/image'
sp.plot.plot_pattern(vmax=99,
                     heatmap=False, 
                     s=5, 
                     reverse_y=True, # optional
                     reverse_x=True, # optional
                     image_path=img_path, # optional
                     rotate_img=True, # optional
                     k=4, # optional
                     aspect=0.55 # optional
                     )
```

<div  align="center">    
  <img src="./pic/scatterplot.png" width = "400" align=center />
</div>

### Visualize the intersections between patterns 3 & 1:

```python
sp.plot.plot_intersection(pattern_list=[0, 1],
                          image_path=img_path,
                          reverse_y=True,
                          reverse_x=True,
                          aspect=0.55,
                          s=20)
```

<div  align="center">    
  <img src="./pic/scatterplot_mx.png" width = "200" align=center />
</div>

### To visualize the gene expression by labels:

```python
sp.plot.plot_genes(label=0, vmax=99)
```

## Attributes of STMiner.SPFinder Object

| Attributes           | Type          | Description                            |
| -------------------- | ------------- | -------------------------------------- |
| adata                | Anndata       | Anndata for loaded spatial data        |
| patterns             | dict          | Spatial distributions pattern of genes |
| genes_patterns       | dict          | GMM model for each gene                |
| global_distance      | pd. DataFrame | Distances between genes and background |
| mds_features         | array         | embedding features of genes            |
| genes_distance_array | pd. DataFrame | Distance between each GMM              |
| genes_labels         | pd. DataFrame | Gene name and their pattern labels     |
| plot                 | Object        | Call plot to visualization             |

# 📜 Release history
https://pypi.org/project/STMiner/#history

# 🔖 Referance
[1] 

# ✉️ Contact
 - Peisen Sun: 📧(sunpeisen@stu.xjtu.edu.cn) / 𝕏(https://x.com/Sun_python)
 - Kai Ye: 📧(kaiye@xjtu.edu.cn)

<br>
Please ⭐Star STMiner on Github if you find it's useful, thank you!
  
<img src="./pic/footer.svg" width="100%">
