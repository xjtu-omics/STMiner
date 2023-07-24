# %%

import numpy as np

a = np.ones((100, 100))

a
# %%
from testUtils import *

b = add_salt_noise(a, 5)
b
# %%

b = np.random.random(a.shape)
import matplotlib.pyplot as plt
import seaborn as sns

sns.distplot(b)
plt.show()

#%%

c = np.random.normal(loc=2,scale=1 ,size=a.shape)
sns.distplot(c)
plt.show()

#%%
c[c<0] = 0
sns.distplot(c)
plt.show()