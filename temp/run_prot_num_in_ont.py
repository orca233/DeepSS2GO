'''
README：用于计算table1，比如 swissprot, ALL的train或test里，计算包含bp/cc/mf的蛋白质分别有多少条

'''

import pandas as pd
import numpy as np

fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/temp/'

df_bp = pd.read_pickle(fpath + 'terms_bp.pkl')
df_cc = pd.read_pickle(fpath + 'terms_cc.pkl')
df_mf = pd.read_pickle(fpath + 'terms_mf.pkl')

df_train = pd.read_pickle(fpath + 'train_data.pkl')
df_test = pd.read_pickle(fpath + 'test_data.pkl')


# 修改下面两行， 分别交叉计算下面 3*2 的结果
df_terms = df_mf
df_dataset = df_train



k = 0
prot_cross = 0
for index, row in df_dataset.iterrows():
    k += 1
    print(k)
    prop_anno = row['prop_annotations']  #
    ont_x = set(df_terms['terms'])  # terms变set
    # print('lol = ', prop_anno & ont_x)

    if len(prop_anno & ont_x) != 0:  # 取交集
        prot_cross += 1
        print(row['proteins'], ', the protein contain this ont')
    else:
        print(row['proteins'], ', the protein NOT contains this ont !!!!!!!!!!!')


print('lol')
print('total protein num = ', k)
print('containing this ont protein num = ', prot_cross)


print('hhh')
