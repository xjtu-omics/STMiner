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

# %%

c = np.random.normal(loc=2, scale=.5, size=a.shape)
sns.distplot(c)
plt.show()

# %%
c[c < 0] = 0
sns.distplot(c)
plt.show()
#%%

a = np.array([[1,2,3],[4,5,6]])
b = np.array([[1,0,1],[0,1,1]])
a*b