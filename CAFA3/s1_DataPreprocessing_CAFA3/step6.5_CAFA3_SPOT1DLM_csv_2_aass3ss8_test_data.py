
# README:
# /SPOT1DLM_results/ 70,000个csv文件（每一个是一个蛋白，包含AA/SS3/SS8/溶解度等信息）提取 AA/SS3/SS8信息，
# 生成新的pkl文件，SPOT1DLM_prot_AA_SS3_SS8.pkl，类似于
# 方便下一步swissprot的替换

import pandas as pd
import os
from step0_DataPreprocessingSetting import *

df_output = pd.DataFrame(columns=['proteins', 'AA', 'SS3', 'SS8'])  # 输出文件4列
proteins = []
AA = []
SS3 = []
SS8 = []

fpath = path_redundancy + "SPOT1DLM_results_CAFA3_test_data/"


fnames = os.listdir(fpath)  # 这是一个list，输出： ['ATD3A_DROME.csv', 'GMC2_YEAST.csv', 'YFFH_SCHPO.csv'...]
#fnames = ['ATD3A_DROME.csv', 'GMC2_YEAST.csv', 'YFFH_SCHPO.csv']  # 仅做尝试，运行时候注释掉！！！！

print('frame length=' + str(len(fnames)))

# 每个文件打开，分别把 名字，AA,SS3,SS8 附加给每个list，最后把list叠加给df_output
k = 0
for fname in fnames:  # fname = file name
    print(k)
    k += 1
    # print(fname)

    prot_df = pd.read_csv(fpath + fname)
    prot_name = fname.split(sep='.')
    proteins.append(prot_name[0])  # name='YFFH_SCHPO.csv', name[0]为前缀

    AA.append(''.join(prot_df['AA'].tolist()))
    SS3.append(''.join(prot_df['SS3'].tolist()))
    SS8.append(''.join(prot_df['SS8'].tolist()))  # prot_df['SS8'].tolist() = ['H', 'H', 'H', 'H', 'H', 'H', 'H', ...]


df_output['proteins'] = proteins
df_output['AA'] = AA
df_output['SS3'] = SS3
df_output['SS8'] = SS8

print('df_output')
print(df_output)

df_output.to_pickle(path_pub_data + 'data_cafa3/CAFA3_test_data_SPOT1DLM_aass3ss8.pkl')  # 包含4列：prot_name, AA, SS3, SS8， 一共69252行，每行一个蛋白质
# 修改：data_cafa3/CAFA3_train_data_SPOT1DLM_aass3ss8.pkl！！！！！！！！！！


# SPOT1DLM_AASS3SS8.pkl一共有四列: 格式
# 第一行是标题： proteins, AA, SS3, SS8
# 之后一共70,000行，每行一个蛋白质
# 'YFFH_SCHPO'
#  'MAILPEFISQTPPVTRYIVLGTLFTTLAVNFGYVSDLKIFFNWKLFLAKGEYWRAITTFLYVGPFGLELILYLSFLLRFMSMLERSSPPPQTQSFLKTVLIVWFSLLVTSYFSYMPFAASYFSFTMLYIWSWKHPLYRISILGLFDVKAPYVPWVMVLLRWLRTGIFPLLDLISALIGHVYFFVTDFSTV',
# 'CCCHHHHHHCCCHHHHHHHHHHHHHHHHHHCCCCCHHHHCECHHHHHCCCCHHHHHHCCCCCCCCCHHHHHHHHHHHHHHHHHHHCCCCCCHHHHHHHHHHHHHHHHHHHHHCCCCCCCHHHHHHHHHHHHHHCCCCEEEECCCEEEEHHHHHHHHHHHHHHHHCCCCHHHHHHHHHHHHHHHHHCCCCC',
# 'CCCHHHHHHTSCHHHHHHHHHHHHHHHHHHTTSSCHHHHCCCHHHHHTTCCHHHHHHTTTCCSSCCHHHHHHHHHHHHHHHHHHHHSSTTCHHHHHHHHHHHHHHHHHHHHHSCCCCCCHHHHHHHHHHHHHHCTTCEEEETTTEEEEGGGHHHHHHHHHHHHHSSCCHHHHHHHHHHHHHHHHHTTCCC']



print('step6_done')




