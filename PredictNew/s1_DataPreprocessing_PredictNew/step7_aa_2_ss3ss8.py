import pandas as pd
from step0_DataPreprocessingSetting import *

swissprot_clean_ALL00_aa_df = pd.read_pickle(path_pub_data + 'data_new/new_clean_aa.pkl')
SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'data_new/new_SPOT1DLM_aass3ss8.pkl')
print(swissprot_clean_ALL00_aa_df['proteins'])
print(SPOT1DLM_aass3ss8)
SPOT1DLM_aass3ss8_dict = SPOT1DLM_aass3ss8.set_index('proteins').T.to_dict('list')

for aa_ss in ['ss8', 'ss3']:
    k = 0
    for index, row in swissprot_clean_ALL00_aa_df.iterrows():
        print(k)
        k += 1
        ss_replace = []
        if aa_ss == 'ss3':
            ss_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][1]  # [0]:AA, [1]:SS3, [2]:SS8
        elif aa_ss == 'ss8':
            ss_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][2]  # [0]:AA, [1]:SS3, [2]:SS8
        swissprot_clean_ALL00_aa_df.loc[index, 'sequences'] = ss_replace

    temp_df = swissprot_clean_ALL00_aa_df

    print('------------ check ------------ ')
    print('new_x_SPOT1DLM_all_ssX_df :::')
    print(temp_df['sequences'])

    temp_df.to_pickle(path_pub_data + 'data_new/new_clean_' + aa_ss + '.pkl')

print('done')
