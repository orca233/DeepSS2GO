"""
记得改step0_setting 中的local setting, train_data = 'all'/MOUSE/ECOLI..... change !!!!!!!!!!!!!!!!!!!!!!!!!!!!

本py改编自：DeepGOPlus中原名：deepgoplus_data.py

第一步：找训练集和测试集
根据train_data (@ step0_TrainTestSetting_local) 是否为ALL00/HUMAN/ECOLI...，看是case 0 or 1 or 3(CAFA3)
case 0:
if train_data == all:  # 训练集=测试集，所以95%的数据用于train, 5%用于test
先有的all_data_df (swissprot_clean_ALL00_aa/ss.pkl)，后分裂出train_data & test_data
其中train_data下一步step11_TrainTest会被split成train & validation

case 1:
elif train_data != test_data: case 1:
先有的train_data & test_data，后合并成all_data_df
cp MOUSE/ECOLI_aa/ss.pkl 得到 train_data/test_data，

case 2:
elif train_data = test_data = CAFA3


第二步：找 terms，满足两点
1 大于 > gominre 默认最少出现50次，得到 terms_GOMinRepeat
2 train & test 中 prop_annotation 交集， 得到 terms_trainXtest
3 上述两个的交集 = terms_GOMinRepeat_trainXtest

P.S. 只满足条件1的terms，拎出来，用于CAFA3_round4，非常还原deepgoplus的套路

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
terms_gominre.pkl  (超过出现50次的， 只满足条件1)
terms_gominre_trxte.pkl  (超过出现50次的 & 同时在train_data.pkl & test_data.pkl 中出现，同时满足条件 1 & 2)

"""

##!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from collections import Counter
# from utils import Ontology, FUNC_DICT
import logging
import os
# import sys
from utils import Ontology, NAMESPACES
from step0_TrainTestSetting_local import *
from step0_TrainTestSetting_global import path_base

logging.basicConfig(level=logging.INFO)


os.system('mkdir -p data')
os.system('mkdir -p results')

@ck.command()
# @ck.option('--go-file', '-gf', default='../../pub_data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')  # output
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')  # output

@ck.option('--terms-gominrepeat-file', '-tgf', default='data/terms_gominre.pkl', help='XX')  # 针对CAFA3_round4，只满足条件1 !!!!
@ck.option('--terms-gominrepeat-trainxtest-file', '-tgtf', default='data/terms_gominre_trxte.pkl', help='XX')  # 同时满足条件1 & 2
@ck.option('--terms-bp-file', '-tbf', default='data/terms_bp.pkl', help='XX')  # 从 terms_GOMinRepeat_trainXtest_file 分化出来（同时满足条件1&2）
@ck.option('--terms-cc-file', '-tcf', default='data/terms_cc.pkl', help='XX')
@ck.option('--terms-mf-file', '-tmf', default='data/terms_mf.pkl', help='XX')

@ck.option('--train-data', '-trd', default=params_local['train_data'], help='HUMAN/ECOLI/...')
@ck.option('--test-data', '-ted', default=params_local['test_data'], help='HUMAN/ECOLI/...')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--gominrepeat', '-gmr', default=params_local['GOMinRepeat'], help='GO min repeat times')
@ck.option('--path-base', '-pb', default=params_local['path_base'], help='..')


def main(go_file, train_data_file, test_data_file, terms_gominrepeat_file, terms_gominrepeat_trainxtest_file,
         terms_bp_file, terms_cc_file, terms_mf_file, train_data, test_data, aa_ss, gominrepeat, path_base):
    print('\n################## a long, long time ago ... ##################\n')
    print('# starting step1_SplitTrainTest_Terms #')
    print('# params_local:')
    for key, value in params_local.items():
        print(key, ' = ', value)

    #########################################
    ######### pre-setting ###################
    #########################################
    path_pub_data = path_base + 'pub_data/'
    path_redundancy = path_base + 'redundancy/'
    # path_s2_TrainTest = path_base + 's2_TrainTest/'


    #####################################################
    ###### step1.1 找训练集和测试集 train & test file ######
    #####################################################

    print('\n-------- Step 1.1: Processing train_data_file & test_data_file ------------')  # logging.info
    # 下面分三类，
    # 1）训练集 = 测试集 != CAFA3
    # 2) 训练集 != 测试集
    # 3） 训练集 = 测试集 = CAFA3

    # 训练集 = 测试集，且 != CAFA3, 比如 TrainHUAMN_TestHUMAN，或者TrainALL00_TestALL00，此时要将所有数据拆成 95% (train+valid) + 5% (test)
    if train_data == test_data and train_data != 'CAFA3':  # train_data = test_data = ALL00/HUMAN/...
        print('\n----- case 0 :train_data == test_data != CAFA3 --------')  # logging.info
        # 数据来源： swissprot_clean_ALL00_aa.pkl / swissprot_clean_HUMAN_ss8.pkl
        swissprot_data_file = path_pub_data + 'swissprot_clean_%s_%s.pkl' % (train_data, aa_ss)
        all_data_df = pd.read_pickle(swissprot_data_file)

        ## spliting all_data to train_data & test_data.pkl
        # 先有的all_data_df，后分裂出train_data & test_data

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


    elif train_data != test_data and train_data != 'CAFA3' and test_data != 'CAFA3':  # train和test data 不一样。比如TrainHUMAN_TestECOLI
        ############## cp special species as train/test_data.pkl ###########
        # 先有的train_data & test_data，后合并成all_data_df
        print('\n----- case 1: train_data !!!=== test_data ------')  # logging.info
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


    elif train_data == test_data and train_data == 'CAFA3':  # 针对CAFA3，train_data & test_data已经确定好了，只需要把terms搞定即可
        print('\n------ case 2: train_data == test_data == CAFA3 -------')  # logging.info
        # 数据来源： CAFA3_train_data_clean_aa/ss8.pkl / CAFA3_test_data_clean_aa/ss8.pkl
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


    ##############################################
    ########## step1.2 找公共 terms ###############
    ##############################################

    print('\n-----Step 1.2: --------- Processing terms ---------------------')  # logging.info

    print('GO Min Repeat=', gominrepeat)  # 这会儿打印出来  有用了！

    ### 先满足条件1： 累积 train + test中，terms出现次数 > gominre=50
    # go = Ontology(go_file, with_rels=True)  # go好像没用上
    cnt = Counter()
    # annotations = list()  # 也没用上？
    for i, row in all_data_df.iterrows():  # original df.iterrows，改成了从训练中找出来的
        for term in row['prop_annotations']:  # 'prop_annotations'是包括从exp_annotation推测出来的
            cnt[term] += 1

    print('Number of terms (raw total = train + test) = ', len(cnt))

    # Filter terms with annotations more than gominre
    res = {}
    for key, val in cnt.items():  # key = 'GO:0072060'
        if val >= gominrepeat:
            ont = key.split(':')[0]  # key.split(':')[0] = 'GO'
            if ont not in res:  # 排除调命名不以GO开头的？ 或者是给res字典开个头？？？
                res[ont] = []
            res[ont].append(key)
    # res 最后是一个字典{'GO': ['GO:0043656', 'GO:0030430', 'GO:0033643', ... ]} 而这里只有一个key，也就是GO

    terms = []  # 筛选后的terms
    for key, val in res.items():
        # print(key, len(val))  # 没有改，因只有GO在terms，这行打出来的就是GO 5242，和下面‘umber of terms (> gominre)’是一样的
        terms += val

    terms_GOMinRepeat = terms  # 这是满足条建1，出现次数的大于 GOMinRe 的terms
    print(f'Number of terms (> gominre) = : {len(terms_GOMinRepeat)}')


    # P.S. CAFA3_round4中terms只考虑terms，不考虑是否和test做交集。i.e. 只考虑条件1，不考虑条件2：（非常还原deepgoplus的套路）

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
    terms_trainXtest = train_x_test_prop_annotations
    print('Number of terms (trainXtest) =', len(terms_trainXtest))

    ### 条件 1 & 2 交集
    terms_GOMinRe_trXte = terms_trainXtest & set(terms_GOMinRepeat)  # 条件1 & 条件2 交集，2088=(1517+269+302)
    print('Number of terms ( >gominre & train_x_test_terms ) = ', len(terms_GOMinRe_trXte))  # = terms_all



    # Save the list of terms
    # 只满足条件1，> GOminrepeat50次  （专门为 CAFA3_round4 准备，非常还原deepgoplus的套路）
    terms_GOMinRepeat = set(terms_GOMinRepeat)
    terms_df_GOMinRepeat = pd.DataFrame({'terms': list(terms_GOMinRepeat)})
    terms_df_GOMinRepeat.to_pickle(terms_gominrepeat_file)

    # 同时满足条件1 & 2，更general,为一般情况准备
    terms_df_gominre_trxte = pd.DataFrame({'terms': list(terms_GOMinRe_trXte)})
    terms_df_gominre_trxte.to_pickle(terms_gominrepeat_trainxtest_file)  # terms 是从’prop_annot‘ 筛出来，重复次数 > GOMinRepeat && 同时出现在train/test_data.pkl中


    ########## 保存 terms_bp/cc/mf.pkl #######################
    # bp=30039, cc=4293, mf=12407
    go_rels = Ontology(go_file, with_rels=True)

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['bp'])  # 所有在当前属于：bp/cc/mf下的GO？
    terms_bp_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_bp_df.to_pickle(terms_bp_file)
    print('terms_bp len=', len(terms_bp_df))  # 1517

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['cc'])
    terms_cc_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_cc_df.to_pickle(terms_cc_file)
    print('terms_cc len=', len(terms_cc_df))  # 269

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['mf'])
    terms_mf_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_mf_df.to_pickle(terms_mf_file)
    print('terms_mf len=', len(terms_mf_df))  # 302



    #  得想办法重新mix一下，最好把 s1 中的 aa/ss3/ss8 合并起来 ！！！！！！！！！！！！！！！！！！！
    #  得想办法重新mix一下，最好把 s1 中的 aa/ss3/ss8 合并起来 ！！！！！！！！！！！！！！！！！！！
    #  得想办法重新mix一下，最好把 s1 中的 aa/ss3/ss8 合并起来 ！！！！！！！！！！！！！！！！！！！

    # 想办法，把train & valid & test 中的terms全部 shuffle 一下！！！！！
    # 想办法，把train & valid & test 中的terms全部 shuffle 一下！！！！！
    # 想办法，把train & valid & test 中的terms全部 shuffle 一下！！！！！


    print('\n################## And they all lived happily ever after! ##################\n')

if __name__ == '__main__':
    main()






