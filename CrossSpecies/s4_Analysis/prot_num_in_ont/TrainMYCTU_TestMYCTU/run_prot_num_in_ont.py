'''
README：用于计算table1，比如 swissprot, ALL的train或test里，计算包含bp/cc/mf的蛋白质分别有多少条

CAFA3/ALL/HUMAN/ECOLI/....每一个类别在一个单独文件夹，该文件夹下的terms_ont & train/test_data.pkl是从相应output/*/data/文件夹cp过来

'''

import pandas as pd
import numpy as np

fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/CrossSpecies/s4_Analysis/prot_num_in_ont/TrainMYCTU_TestMYCTU/'  # 修改！！！！！！

df_bp = pd.read_pickle(fpath + 'terms_bp.pkl')
df_cc = pd.read_pickle(fpath + 'terms_cc.pkl')
df_mf = pd.read_pickle(fpath + 'terms_mf.pkl')

df_train = pd.read_pickle(fpath + 'train_data.pkl')
df_test = pd.read_pickle(fpath + 'test_data.pkl')

# 修改下面两行， 分别交叉计算下面 3*2 的结果  # 修改！！！！！！
df_terms = df_mf  # df_bp, df_cc, df_mf
df_dataset = df_test  # df_train, df_test

k = 0
prot_with_ont = 0
for index, row in df_dataset.iterrows():
    k += 1
    print(k)
    prop_anno = row['prop_annotations']  #
    if isinstance(prop_anno, list):
        prop_anno = set(prop_anno)
    ont_x = set(df_terms['terms'])  # terms变set
    # print('lol = ', prop_anno & ont_x)

    if len(prop_anno & ont_x) != 0:  # 取交集
        prot_with_ont += 1
        print(row['proteins'], ', the protein contain this ont')
    else:
        print(row['proteins'], ', the protein NOT contains this ont !!!!!!!!!!!')


print('lol')
print('protein num total = ', k)
print('protein num containing this ont = ', prot_with_ont)
print('terms classes = ', len(df_terms['terms']))




