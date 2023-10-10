import pandas as pd
import numpy as np

kernels = [[16, 16], [32, 32], [32, 32]]
filters = [[65536, 65536]]

for i in kernels:
    for j in filters:
        print(i, 'lol', j)

x = np.array(kernels)
print(x.shape)
print(type(x.shape))

y = np.array([16, 32, 64])
print(y.shape)
print(type(y.shape))
