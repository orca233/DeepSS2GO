
### 重要 !!!!!!!!!!!!!!!
# README
# 从DeepGOPlus step1 得到的swissprot.pkl  因为SPOT1DLM的两个部分，esm和prottrans对序列有要求，所以进行清理数据
# seq长度超过1024-esm不开心/
# 超过900-GPUcuda内存不开心/ ----这个在v100可以不考虑？？？
# seq和aa并集长度!=20,有其他元素，不开心

# OUT： swissprot_clean_ALL00_aa.pkl  经过SPOT1DLM优化清洗过的swissprot
# 这里生成的是ALL所有prot的seq aa序列，下一步step345再通过SPOT1DLM预测出对应的ss

import pandas as pd
from step0_DataPreprocessingSetting import *

df = pd.read_pickle(path_pub_data + 'swissprot.pkl')  # 全集 79230
print('original swissprot.pkl::::::::::::::::::')
print(df.info())
print(df)

aa = {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y'}

k = 0
for index, row in df.iterrows():
    k += 1
    print(k)
    seq = row['sequences']
    if len(seq) > 1021 or len(aa | set(seq)) != 20:  # seq长度超过1024/ 超过900- 或者 seq和aa并集长度!=20,有其他元素，则删除
        df = df.drop(index=index)  # 删除命令

# output
df.to_pickle(path_pub_data + 'swissprot_clean_ALL00_aa.pkl')  # final  pkl = 79,020

print('swissprot_clean_ALL00_aa.pkl:::::')
print(df.info())
print(df)


print('step2 done')




