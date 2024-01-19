
import pandas as pd
from step0_DataPreprocessingSetting import *


df = pd.read_pickle(path_base + 'pub_data/data_new/new_aa.pkl')
print('original new_aa.pkl::::::::::::::::::')
print(df.info())
print(df)

aa = {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y'}

k = 0
for index, row in df.iterrows():
    k += 1
    print(k)
    seq = row['sequences']
    if len(seq) > 1021 or len(aa | set(seq)) != 20:
        df = df.drop(index=index)

# output
df.to_pickle(path_base + 'pub_data/data_new/new_clean_aa.pkl')

print('new_clean_aa.pkl:::::')
print(df.info())
print(df)


print('step2 done')




