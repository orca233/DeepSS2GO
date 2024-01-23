import pandas as pd
from step0_DataPreprocessingSetting import *

print('\n################## a long, long time ago ... ##################\n')
print('# starting step2_Swissprot_x_SPOT1DLM #')

df = pd.read_pickle(path_base + 'pub_data/swissprot.pkl')
print('original swissprot.pkl::::::::::::::::::')
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
df.to_pickle(path_pub_data + 'swissprot_clean_ALL00_aa.pkl')

print('swissprot_clean_ALL00_aa.pkl:::::')
print(df.info())
print(df)

print('step2 done')

print('\n################## And they all lived happily ever after! ##################\n')

