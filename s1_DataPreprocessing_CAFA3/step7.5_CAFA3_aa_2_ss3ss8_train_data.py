import pandas as pd
# import tensorflow as tf
import os
from step0_DataPreprocessingSetting import *

# README:
# INPUT：
# “swissprot_clean_ALL00_aa.pkl”     (包含10列，prot_name, seq, GO, CAFA等信息，step2_swissprot_x_SPOT1DLM_aa.py生成)
# “SPOT1DLM_aass3ss8.pkl”        (包含4列，prot_name, AA, SS3, SS8对应序列顺序，step6_SPOT1DLM_csv_2_AASS3SS8.py生成)
# OUT:
# 后者的ss3/ss8信息替换给前者的’seq‘，生成'swissprot_clean_ALL00_aa_ss8/ss3.pkl'两个文件


swissprot_clean_ALL00_aa_df = pd.read_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_clean_aa.pkl')  # 另一个: data_cafa3/CAFA3_train_data_clean_aa.pkl
SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_SPOT1DLM_aass3ss8.pkl')  # original: SPOT1DLM_aass3ss8.pkl

print(swissprot_clean_ALL00_aa_df['proteins'])
print(SPOT1DLM_aass3ss8)

# swissprot_clean_ALL00_aa_df = swissprot_clean_ALL00_aa_df[0:100] # 测试前100行，正式运行时注释掉

# df转成以一列为key的dict
SPOT1DLM_aass3ss8_dict = SPOT1DLM_aass3ss8.set_index('proteins').T.to_dict('list')


# SPOT1DLM_AASS3SS8_dict: 格式
# {第一个  ... 'GMC2_YEAST': ['AA', 'SS3', 'SS8'], 'YFFH_SCHPO' : ['AA', 'SS3', 'SS8']}

# 'YFFH_SCHPO': [
#  'MAILPEFISQTPPVTRYIVLGTLFTTLAVNFGYVSDLKIFFNWKLFLAKGEYWRAITTFLYVGPFGLELILYLSFLLRFMSMLERSSPPPQTQSFLKTVLIVWFSLLVTSYFSYMPFAASYFSFTMLYIWSWKHPLYRISILGLFDVKAPYVPWVMVLLRWLRTGIFPLLDLISALIGHVYFFVTDFSTV',
# 'CCCHHHHHHCCCHHHHHHHHHHHHHHHHHHCCCCCHHHHCECHHHHHCCCCHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHHHHHHHHHHCCCCCCCHHHHHHHHHHHHHHCCCCEEEECCCEEEEHHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHCCCCC',
# 'CCCHHHHHHTSCHHHHHHHHHHHHHHHHHHTTSSCHHHHCCCHHHHHTTCCHHHHHHTTTCCSSCCHHHHHHHHHHHHHHHHHHHHSSTTCHHHHHHHHHHHHHHHHHHHHHSCCCCCCHHHHHHHHHHHHHHCTTCEEEETTTEEEEGGGHHHHHHHHHHHHHSSCCHHHHHHHHHHHHHHHHHTTCCC']


# 逐行确定swissprot_clean_aa_df中的row['proteins'] -- prot_name，反定位在SPOT1DLM_AASS3SS8_dict中找到SS8，替换’seq‘

for aa_ss in ['ss8', 'ss3']:

# aa_ss = 'ss3'  # change!!!!!!!!!!!!!!!!!!!!!!!! ss3 or ss8

    k = 0  # 计数
    for index, row in swissprot_clean_ALL00_aa_df.iterrows():
        print(k)
        k += 1
        ss_replace = []  # 这个ss是要用于替换 swissprot_x_SPOT1DLM_all_aa 中 seq 的
        if aa_ss == 'ss3':
            ss_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][1]  # [0]:AA, [1]:SS3, [2]:SS8
        elif aa_ss == 'ss8':
            ss_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][2]  # [0]:AA, [1]:SS3, [2]:SS8
        swissprot_clean_ALL00_aa_df.loc[index, 'sequences'] = ss_replace  # 定位&替换成test_ss中的对应第i个二级结构序列  # 替换


    temp_df = swissprot_clean_ALL00_aa_df   # 这个可能是ss3或ss8替换了aa之后的df

    print('------------ check ------------ ')
    print('swissprot_x_SPOT1DLM_all_ssX_df :::')
    print(temp_df['sequences'])  # 看看替换了ss3，ss8没有

    # 存成swissprot_clean_ALL00_aa_ss8.pk  &  swissprot_clean_ALL00_aa_ss3.pk
    temp_df.to_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_clean_' + aa_ss + '.pkl')  # 第二遍用： CAFA3_clean_test_' + aa_ss + '.pkl



print('done')
