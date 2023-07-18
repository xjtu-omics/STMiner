#%%
from SPFinder import SPFinder

#%%
spf = SPFinder()
spf.read_10x(file='F://Rep11_MOB_ST.h5ad', amplification=1000, bin_size=80)

#%%
spf.normalize()
spf.fit_pattern(100, n_comp=10)
spf.cluster(6)

#%%
spf.genes_labels
