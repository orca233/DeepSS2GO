
### 重要 !!!!!!!!!!!!!!!
# README
# 从DeepGOPlus step1 得到的CAFA3_test_data_raw.pkl & CAFA3_train_data_raw.pkl
# 因为SPOT1DLM的两个部分，esm和prottrans对序列有要求，所以进行清理数据

# seq长度超过1024-esm不开心/
# 超过900-GPUcuda内存不开心/ ----这个在v100可以不考虑？？？
# seq和aa并集长度!=20,有其他元素，不开心

# input:
# CAFA3_train_data_raw.pkl  (56M, 66,841)
# CAFA3_test_data_raw.pkl (2M, 3,328)

# OUT：
# CAFA3_clean_train_data_aa.pkl
# CAFA3_clean_test_data_aa.pkl
# 经过SPOT1DLM优化清洗过的
# 这里生成的是ALL所有prot的seq aa序列，下一步step345再通过SPOT1DLM预测出对应的ss

import pandas as pd
from step0_DataPreprocessingSetting import *



def cut_long_aa(df):
    aa = {'A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y'}
    k = 0
    for index, row in df.iterrows():
        k += 1
        print(k)
        seq = row['sequences']
        if len(seq) > 1021 or len(aa | set(seq)) != 20:  # seq长度超过1024/ 超过900- 或者 seq和aa并集长度!=20,有其他元素，则删除
            df = df.drop(index=index)  # 删除命令

    return df

### train ###
df_train = pd.read_pickle(path_pub_data + 'data_cafa3/CAFA3_train_data_raw.pkl')  #
df_train_clean = cut_long_aa(df_train)
df_train_clean.to_pickle(path_pub_data + 'data_cafa3/CAFA3_train_data_clean_aa.pkl')  #

### test ###
df_test = pd.read_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_raw.pkl')  #
df_test_clean = cut_long_aa(df_test)
df_test_clean.to_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_clean_aa.pkl')  #


print('train_raw\n', df_train.info())  # 66841
print('train_clean\n', df_train_clean.info())  # 60372
print('test_raw\n', df_test.info())  # 3328
print('test_clean\n', df_test_clean.info())  # 3049


print('step2 done')




