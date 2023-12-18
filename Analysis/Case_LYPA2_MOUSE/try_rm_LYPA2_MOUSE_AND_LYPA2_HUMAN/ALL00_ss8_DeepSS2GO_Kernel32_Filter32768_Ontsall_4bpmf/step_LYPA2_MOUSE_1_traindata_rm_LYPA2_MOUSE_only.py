
# 首先把/data/中的 train_data.pkl & test_data.pkl 备份为 train_data_bcp.pkl & test_data_bcp.pkl
#
# 查看在 data/中的 train_data.pkl & test_data.pkl是否含有LYPA2_MOUSE这个蛋白。如果存在，删除。
# 更新的 train_data.pkl，删除LYPA2_MOUSE & LYPA2_HUMAN



import pandas as pd
import numpy as np


# dir_path = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/Analysis/Case_LYPA2_MOUSE/'
#
# fpath_train = dir_path + 'ALL00_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall_4bpccmf/data/train_data_bcp.pkl'
# fpath_test = dir_path + 'ALL00_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall_4bpccmf/data/test_data_bcp.pkl'


fpath_train = './data/train_data_bcp.pkl'
fpath_test = './data/test_data_bcp.pkl'


df_train = pd.read_pickle(fpath_train)
df_test = pd.read_pickle(fpath_test)

print(df_train.info())



# 检查 train_data_bcp.pkl 是否包含 LYPA2： 包含的，一个HUMAN，一个MOUSE
print('---- check whether train/test_data.pkl contain LYPA2 ----')

for i in df_train['proteins']:
    if 'LYPA2' in i:
        print(i)


# 删除
print('----- recheck and del -----')
k = 0
for index, row in df_train.iterrows():
    k += 1
    # print(k)
    if 'LYPA2_MOUSE' in row['proteins']:  # 不同于另一个，会删除所有的LYPA2_MOUSE&HUMAN
        print(row['proteins'])
        df_train = df_train.drop(index=index)  # 删除命令

print('k=', k)

# output
df_train.to_pickle('./data/train_data.pkl')  # final  pkl = 79,020

# 看看新的train_data.pkl有多少条蛋白质（旧的有68325条蛋白质）

df_new = pd.read_pickle('./data/train_data.pkl')
print('updated new train_data contains proteins number = ?')  # 68323
print(df_new.info())




