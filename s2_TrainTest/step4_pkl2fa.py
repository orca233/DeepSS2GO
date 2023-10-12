"""
README
分几个步骤

第一步，为Diamond做数据预处理，得到 train/test_data.fa文件
case 0) if CAFA3
case 1) elif aa_ss = aa, 分三步：pkl2fa, build *dmnd, querying 生成 diamond_aa.res
case 2) elif aa_ss = ss8/ss3，则从 *aa/data/diamond_aa.res cp 到 *ss/data/.

# 第二步，Diamond  (不论 aa / ss，都是用 aa 做diamond比对)
#
# -d          = train_data.dmnd (数据库)
# -q (query)  = test_data.fa （要查询的文件）
# -o (output) = test_diamond.res （输出文件）
# qseqid      = Query Seq-id    查询列
# qseqid      = Subject Seq-id  数据库列
# bitscore    = Bit score     分数
#
# IN:
# train_data.dmnd    Diamond数据库
# test_data.fa            查询文件
#
# OUT:
# diamond_aa.res
"""


import click as ck
import os
import pandas as pd
# from step0_TrainTestSetting_global import *
from step0_TrainTestSetting_local import *  # 里面有aa_ss信息
from utils import pkl2fa

@ck.command()
# @ck.option('--go-file', '-gf', default='../../pub_data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')

@ck.option('--terms-gominrepeat-file', '-tgf', default='data/terms_gominre.pkl', help='XX')  # 针对CAFA3_round4，只满足条件1
@ck.option('--terms-gominrepeat-trainxtest-file', '-tgtf', default='data/terms_gominre_trxte.pkl', help='XX')  # 同时满足条件1 & 2
@ck.option('--terms-bp-file', '-tbf', default='data/terms_bp.pkl', help='XX')  # 从 terms_GOMinRepeat_trainXtest_file 分化出来（同时满足条件1&2）
@ck.option('--terms-cc-file', '-tcf', default='data/terms_cc.pkl', help='XX')
@ck.option('--terms-mf-file', '-tmf', default='data/terms_mf.pkl', help='XX')

@ck.option('--train-data', '-trd', default=params_local['train_data'], help='HUMAN/ECOLI/...')
@ck.option('--test-data', '-ted', default=params_local['test_data'], help='HUMAN/ECOLI/...')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--gominrepeat', '-gmr', default=params_local['GOMinRepeat'], help='GO min repeat times')
@ck.option('--path-base', '-pb', default=params_local['path_base'], help='..')


def main(train_data_file, test_data_file, terms_gominrepeat_file, terms_gominrepeat_trainxtest_file,
         terms_bp_file, terms_cc_file, terms_mf_file, train_data, test_data, aa_ss, gominrepeat, path_base):

    print('\n################## a long, long time ago ... ##################\n')
    print('# starting step4_FindAlpha #')
    print('\n-------- Step 1.1: Processing train_data.fa & test_data.fa ------------')  # logging.info
    print('filename_path = ', os.getcwd())

    path_pub_data = path_base + 'pub_data/'

    ################################################################################
    ############ 第一步，为diamond做数据预处理，得到 train/test_data.fa 文件 #############
    ################################################################################

    # case 0, CAFA3， 直接cp过来， 两个fa文件已经在s1中预处理好了
    if params_local['train_data'] == 'CAFA3' and params_local['test_data'] == 'CAFA3':  # 对于CAFA3，因为train & test 的aa/ss8 已经拿到fa了，所以cp过来直接diamond
        print('case_0: train_data = CAFA3')
        os.system(f'cp {path_pub_data}/data_cafa3/CAFA3_train_data_clean_aa.fa data/train_data.fa')
        os.system(f'cp {path_pub_data}/data_cafa3/CAFA3_test_data_clean_aa.fa data/test_data.fa')

    # case 1, 'aa'，此时需要用 pkl2fa 来将原始的 train_data.pkl 转成 train_data.fa
    elif params_local['aa_ss'] == 'aa':  # params_local['train_data'] != 'CAFA3'
        print('---------- case 1, aa_ss = aa, preparing test_diamond.res ----------')

        # pkl2fa
        print('pkl-2-fa starting ...') # 修改！！！！！！！！！！！！！！
        pkl2fa(train_data_file, 'data/train_data.fa')  # data/train_data.pkl
        pkl2fa(test_data_file, 'data/test_data.fa')  # data/test_data.pkl

        # os.system('python Alpha_pkl2fa.py -df data/train_data.pkl -o data/train_data.fa')  #
        # os.system('python Alpha_pkl2fa.py -df data/test_data.pkl -o data/test_data.fa')

    # case 2, 'ss8'
    elif params_local['aa_ss'] == 'ss8' or params_local['aa_ss'] == 'ss3':  # aa_ss = ss8 or ss3 & params_local['train_data'] != 'CAFA3'
        print('--------- case 2: aa_ss = ss8 or ss3, PLEASE run FS_ss8_2_aa first --------------')

        ################## FS_ss8_2_aa ##########################
        # 将 train_data.pkl & test_data.pkl 转换成 train_data_aa.pkl & test_data_aa.pkl， sequence这一列从ss8/ss3变回aa
        # 原来的是ss8/ss3，现在的sequence一列变成aa。这是为了diamond使用
        train_data_ss8_df = pd.read_pickle(train_data_file)  # original: './data/train_data.pkl'
        test_data_ss8_df = pd.read_pickle(test_data_file)  # orignal: './data/test_data.pkl'
        SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'SPOT1DLM_aass3ss8.pkl')

        print('\ntrain_data_ss8_df = \n', train_data_ss8_df['proteins'])
        print('\ntest_data_ss8_df = \n', test_data_ss8_df['proteins'])
        print('\nSPOT1DLM_aass3ss8 = \n', SPOT1DLM_aass3ss8)

        # df转成以一列为key的dict
        SPOT1DLM_aass3ss8_dict = SPOT1DLM_aass3ss8.set_index('proteins').T.to_dict('list')

        # SPOT1DLM_AASS3SS8_dict: 格式
        # {第一个  ... 'GMC2_YEAST': ['AA', 'SS3', 'SS8'], 'YFFH_SCHPO' : ['AA', 'SS3', 'SS8']}

        def ss8_2_aa(df):
            k = 0  # 计数
            for index, row in df.iterrows():  # swissprot_clean_ALL00_aa_df
                # print(k)
                k += 1
                aa_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][0]  # [0]:AA, [1]:SS3, [2]:SS8
                df.loc[index, 'sequences'] = aa_replace  # 定位&替换成test_ss中的对应第i个二级结构序列  # 替换
            temp_df = df  # 这个可能是ss3或ss8替换了aa之后的df
            return temp_df

        print('222222222222222222222')
        temp_df = ss8_2_aa(train_data_ss8_df)
        temp_df.to_pickle('./data/train_data_aa.pkl')
        print('train_data_aa.pkl = \n', temp_df.info())
        print(temp_df['sequences'])

        temp_df = ss8_2_aa(test_data_ss8_df)
        temp_df.to_pickle('./data/test_data_aa.pkl')
        print('test_data_aa.pkl = \n', temp_df.info())
        print(temp_df['sequences'])
        # os.system('python FS_ss8_2_aa.py')



        # pkl2fa
        print('pkl-2-fa starting ...')
        pkl2fa('data/train_data_aa.pkl', 'data/train_data.fa')  # data/train_data.pkl  # 这里改动了！！！20230630. 另外train_data_file在前面也不对！
        pkl2fa('data/test_data_aa.pkl', 'data/test_data.fa')  # data/test_data.pkl  test_data_file写在前面也不对！！！
        # os.system('python Alpha_pkl2fa.py -df data/train_data_aa.pkl -o data/train_data.fa')  # 这里改动了！！！20230630
        # os.system('python Alpha_pkl2fa.py -df data/test_data_aa.pkl -o data/test_data.fa')

    else:
        print('Holy crap! No specific cases...')



    #
    # ######################################################################
    # ###################### 第二步，diamond #################################
    # ######################################################################
    #
    # # 建库， train_data.fa--> train_data.dmnd 库
    # print('\n-------------- build database train_data.dmnd --------------')
    # os.system('diamond makedb --in data/train_data.fa -d data/train_data')
    # print('\n')
    #
    # # 比对 查询 querying
    # print('-------------- querying from database, creating diamond_aa.res --------------')
    # os.system('diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp '
    #           '-q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res')


    print('\n################## And they all lived happily ever after! ##################\n')

    ############# end #################




if __name__ == '__main__':
    main()



