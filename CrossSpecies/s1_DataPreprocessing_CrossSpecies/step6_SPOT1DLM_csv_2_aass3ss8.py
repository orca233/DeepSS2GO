import pandas as pd
from step0_DataPreprocessingSetting import *

print('\n################## a long, long time ago ... ##################\n')
print('# starting step6_SPOT1DLM_csv_2_aass3ss8 #')


df_output = pd.DataFrame(columns=['proteins', 'AA', 'SS3', 'SS8'])  # 输出文件4列
proteins = []
AA = []
SS3 = []
SS8 = []

fpath = path_redundancy + "SPOT1DLM_results/"

fnames = os.listdir(fpath)

print('frame length=' + str(len(fnames)))

k = 0
for fname in fnames:
    print(k)
    k += 1

    prot_df = pd.read_csv(fpath + fname)
    prot_name = fname.split(sep='.')
    proteins.append(prot_name[0])

    AA.append(''.join(prot_df['AA'].tolist()))
    SS3.append(''.join(prot_df['SS3'].tolist()))
    SS8.append(''.join(prot_df['SS8'].tolist()))


df_output['proteins'] = proteins
df_output['AA'] = AA
df_output['SS3'] = SS3
df_output['SS8'] = SS8

print('df_output')
print(df_output)

df_output.to_pickle(path_pub_data + 'SPOT1DLM_aass3ss8.pkl')



print('step6_done')
print('\n################## And they all lived happily ever after! ##################\n')


