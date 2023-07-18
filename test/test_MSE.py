#%%
from SPFinder import SPFinder

#%%
spf = SPFinder()
spf.read_10x(file='F://Rep11_MOB_ST.h5ad', amplification=1000, bin_size=20)

#%%
spf.