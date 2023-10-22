
"""
记得改step0_setting 中的local setting, train_data = 'all'/MOUSE/ECOLI..... change !!!!!!!!!!!!!!!!!!!!!!!!!!!!

DeepGOPlus中原名：deepgoplus_data.py

第一步：找训练集和测试集
根据train_data (@ step0_setting) 是否为ALL，看是case 0 or 1
if train_data == all: case 0:  # 训练集=测试集，所以95%的数据用于train, 5%用于test
先有的all_data_df (swissprot_clean_ALL00_aa/ss.pkl)，后分裂出train_data & test_data
其中train_data下一步step11_TrainTest会被split成train & validation

else: case 1:
先有的train_data & test_data，后合并成all_data_df
cp MOUSE/ECOLI_aa/ss.pkl 得到 train_data/test_data，


第二步：找 terms，满足两点
1 大于 > gominre 默认最少出现50次
2 train & test 中 prop_annotation 交集
3 上述两个的交集 = final_terms

对于TrainALL_TestALL，terms共cnt=30285，经过>gominre 50 筛选，剩5242种， 这个数字是基于原始的data_all_df的，而非单一train_data
再和train & test_data.pkl公共prop_annotation求交集，得到5164条terms，即出现在step11_TrainTest后的predictions.pkl中的label列数

#################
IN:
path_SPOT1DLM_data + 'swissprot_x_SPOT1DLM_all_aa.pkl'
go.obo

OUT:
/data/
train_data.pkl
test_data.pkl
terms.pkl  (超过出现50次的 & 同时在train_data.pkl & test_data.pkl 中出现)
"""




##!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from collections import Counter
from utils import Ontology, FUNC_DICT
import logging
import os
import sys
from utils import FUNC_DICT, Ontology, NAMESPACES

from step0_TrainTestSetting_local import *

logging.basicConfig(level=logging.INFO)



print('\n\n\n\n starting step1 #####################')

for key, value in params_local.items():
    print(key, ' = ', value)

print('\n')

#### setting #### split后，重命名为train/test_data 并保存在 data/文件夹下
train_data_file = 'data/train_data.pkl'
test_data_file = 'data/test_data.pkl'
# out_terms_file = 'data/terms.pkl'  # 同时出现在train&test中，重复次数超过GOMinRepeat，

terms_all_file = 'data/terms_all.pkl'  # 同时出现在train&test中，重复次数超过GOMinRepeat，
terms_bp_file = 'data/terms_bp.pkl'  # terms_all.pkl中，属于bp部分
terms_cc_file = 'data/terms_cc.pkl'
terms_mf_file = 'data/terms_mf.pkl'


### 从 step0_TrainTestSetting_local.py中导入的部分变量
train_data = params_local['train_data']
test_data = params_local['test_data']
aa_ss = params_local['aa_ss']
GOMinRepeat = params_local['GOMinRepeat']

path_pub_data = params_local['path_pub_data']
go_file = params_local['go_file']



#########################################
##########  找训练集和测试集 ###############
#########################################

# 训练集 = 测试集
if train_data == test_data and train_data != 'CAFA3':  # 所有数据95% train+valid，5% test，train_data = test_data = ALL00/HUMAN/...
    # 数据来源： swissprot_clean_ALL00_aa.pkl / swissprot_clean_HUMAN_ss8.pkl
    swissprot_data_file = path_pub_data + 'swissprot_clean_%s_%s.pkl' % (train_data, aa_ss)
    all_data_df = pd.read_pickle(swissprot_data_file)

    ## spliting all_data to train_data & test_data.pkl
    # 先有的all_data_df，后分裂出train_data & test_data

    logging.info('-----STEP1: Processing train_data & test_data ----- case 0 :train_data == test_data --------')

    n = len(all_data_df)
    # Split train/valid
    index = np.arange(n)
    train_n = int(n * params_local['TrainTestRatio'])  # original=n*0.95  默认traintestratio=0.95
    np.random.seed(seed=0)  # ？
    np.random.shuffle(index)
    train_df = all_data_df.iloc[index[:train_n]]
    test_df = all_data_df.iloc[index[train_n:]]

    # save pkl
    train_df.to_pickle(train_data_file)
    test_df.to_pickle(test_data_file)

    print('train data = test data = %s' % train_data)
    print("Number of all proteins", len(all_data_df))
    print('Number of train proteins', len(train_df))
    print('Number of test proteins', len(test_df))


elif train_data != test_data and train_data != 'CAFA3':
    # train_和 test data 不一样。比如TrainHUMAN_TestECOLI MOUSE/ECOLI...
    ############## cp special species as train/test_data.pkl ###########
    # 先有的train_data & test_data，后合并成all_data_df

    # print('wwww3333333333333333333333')

    logging.info(('-----STEP1: Processing train_data & test_data ----- case 1: train_data !!!=== test_data ---'))  # 这也是step1 wow

    os.system('cp %sswissprot_clean_%s_%s.pkl data/train_data.pkl' % (path_pub_data, train_data, aa_ss))
    os.system('cp %sswissprot_clean_%s_%s.pkl data/test_data.pkl' % (path_pub_data, test_data, aa_ss))

    # read pkl
    train_df = pd.read_pickle('data/train_data.pkl')
    test_df = pd.read_pickle('data/test_data.pkl')
    all_data_df = pd.concat([train_df, test_df]).reset_index(drop=True)

    print('train data = %s' % train_data)
    print('test data = %s' % test_data)
    print("Number of all proteins (containing 2 species)", len(all_data_df))  # 总长
    print('Number of train proteins', len(train_df))
    print('Number of test proteins', len(test_df))


elif train_data == test_data and train_data == 'CAFA3':
    # 针对CAFA3，train_data & test_data已经确定好了，只需要把terms搞定即可

    logging.info(
        ('-----STEP1: Processing train_data & test_data ----- case 2: train_data == test_data == CAFA3 -------'))  # 这也是step1 wow


    os.system('cp %sdata_cafa3/CAFA3_train_data_clean_%s.pkl data/train_data.pkl' % (path_pub_data, aa_ss))
    os.system('cp %sdata_cafa3/CAFA3_test_data_clean_%s.pkl data/test_data.pkl' % (path_pub_data, aa_ss))


    # read pkl
    train_df = pd.read_pickle('data/train_data.pkl')
    test_df = pd.read_pickle('data/test_data.pkl')
    all_data_df = pd.concat([train_df, test_df]).reset_index(drop=True)

    print('train data = %s' % train_data)
    print('test data = %s' % test_data)
    print("Number of all proteins (containing 2 species)", len(all_data_df))  # 总长
    print('Number of train proteins', len(train_df))
    print('Number of test proteins', len(test_df))


#########################################
##########  找公共 terms ###############
#########################################

logging.info('-----SPTE2: --------- Processing terms ---------------------')

print('GO Min Repeat=', GOMinRepeat)  # 这会儿打印出来  有用了！

### 先满足条件1： 大于 > gominre=50 条件

# print('loading GO ...')
go = Ontology(go_file, with_rels=True)  # go好像没用上

# logging.info('Processing annotations ...')

cnt = Counter()
annotations = list()
for i, row in all_data_df.iterrows():  # original df.iterrows，改成了从训练中找出来的
    for term in row['prop_annotations']:  # 'prop_annotations'是包括从exp_annotation推测出来的
        cnt[term] += 1

print('Number of terms (total)=', len(cnt))

# Filter terms with annotations more than gominre
res = {}
for key, val in cnt.items():  # key = 'GO:0072060'
    if val >= GOMinRepeat:
        ont = key.split(':')[0]  # key.split(':')[0] = 'GO'

        if ont not in res:  # 排除调命名不以GO开头的？ 或者是给res字典开个头？？？
            res[ont] = []
        res[ont].append(key)

# res 最后是一个字典{'GO': ['GO:0043656', 'GO:0030430', 'GO:0033643', ... ]} 而这里只有一个key，也就是GO
# print(re-s)

terms = []  # 筛选后的terms
for key, val in res.items():
    # print(key, len(val))  # 没有改，因只有GO在terms，这行打出来的就是GO 5242，和下面‘umber of terms (> gominre)’是一样的
    terms += val

# print(terms)
print(f'Number of terms (> gominre): {len(terms)}')


'''
# 以下部分注释是为了CAFA3_round4中terms只考虑terms，不考虑是否和test做交集。

### 满足条件2：同时在train_data & test_data.pkl中出现
# 分别求出train / test prop_annotations
train_prop_annotation = set()
test_prop_annotation = set()
for index, row in train_df.iterrows():
    prop_annotations = train_df.loc[index, 'prop_annotations']  # 这是个list ['GO:0043228', 'GO:0048856', ...]
    train_prop_annotation = train_prop_annotation | set(prop_annotations)  # 并集 更新
for index, row in test_df.iterrows():
    prop_annotations = test_df.loc[index, 'prop_annotations']  # 这是个list ['GO:0043228', 'GO:0048856', ...]
    test_prop_annotation = test_prop_annotation | set(prop_annotations)  # 并集 更新

train_x_test_prop_annotations = train_prop_annotation & test_prop_annotation  # train & test 注释集 交集
print('train_x_test_terms =', len(train_x_test_prop_annotations))

### 条件 1 & 2 交集
final_terms = train_x_test_prop_annotations & set(terms)  # 上述集合和 > gominre 取交集   2088=(1517+269+302)
print('train_x_test_terms & > gominre -- terms_all -- =', len(final_terms))

'''


# 专门为 CAFA3_round4 准备的下面这一行：
final_terms = set(terms)  # 上述集合和 > gominre 取交集   2088=(1517+269+302)
print('\nfinal_terms =', len(final_terms))





# Save the list of terms
terms_df = pd.DataFrame({'terms': list(final_terms)})
terms_df.to_pickle(terms_all_file)  # terms 是从’prop_annot‘ 筛出来，重复次数 > GOMinRepeat && 同时出现在train/test_data.pkl中

# print(terms_df)
# print(terms_df.info())



########## 保存 terms_all/bp/cc/mf.pkl #######################
# bp=30039, cc=4293, mf=12407
go_rels = Ontology(go_file, with_rels=True)

go_set_ont = go_rels.get_namespace_terms(NAMESPACES['bp'])  # 所有在当前属于：bp/cc/mf下的GO？
terms_bp_df = terms_df[terms_df['terms'].isin(go_set_ont)]
terms_bp_df.to_pickle(terms_bp_file)
print('terms_bp len=', len(terms_bp_df))  # 1517

go_set_ont = go_rels.get_namespace_terms(NAMESPACES['cc'])  
terms_cc_df = terms_df[terms_df['terms'].isin(go_set_ont)]
terms_cc_df.to_pickle(terms_cc_file)
print('terms_cc len=', len(terms_cc_df))  # 269

go_set_ont = go_rels.get_namespace_terms(NAMESPACES['mf'])  
terms_mf_df = terms_df[terms_df['terms'].isin(go_set_ont)]
terms_mf_df.to_pickle(terms_mf_file)
print('terms_mf len=', len(terms_mf_df))  # 302


# [terms_df['terms'].isin(go_set_ont)] 输出布尔值，
# [0 True
#  1 False
#  2 True
#  ...



print('step1 done \n')


# if __name__ == '__main__':
#     main()


