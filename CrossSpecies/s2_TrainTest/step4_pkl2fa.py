import click as ck
import os
import pandas as pd
from step0_TrainTestSetting_local import *
from utils import pkl2fa
from step0_TrainTestSetting_global import path_base


@ck.command()
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')
@ck.option('--path-base', '-pb', default=path_base, help='..')


def main(train_data_file, test_data_file, path_base):
    print('\n################## a long, long time ago ... ##################\n')
    print('# starting step4_FindAlpha #')
    print('\n-------- Step 1.1: Processing train_data.fa & test_data.fa ------------')
    print('filename_path = ', os.getcwd())

    path_pub_data = path_base + 'pub_data/'

    ############ 4 diamond #############

    if params_local['train_data'] == 'CAFA3' and params_local['test_data'] == 'CAFA3':
        print('case_0: train_data = CAFA3')
        os.system(f'cp {path_pub_data}/data_cafa3/CAFA3_train_data_clean_aa.fa data/train_data.fa')
        os.system(f'cp {path_pub_data}/data_cafa3/CAFA3_test_data_clean_aa.fa data/test_data.fa')

    elif params_local['aa_ss'] == 'aa':
        print('---------- case 1, aa_ss = aa, preparing test_diamond.res ----------')

        # pkl2fa
        print('pkl-2-fa starting ...')
        pkl2fa(train_data_file, 'data/train_data.fa')
        pkl2fa(test_data_file, 'data/test_data.fa')

    # case 2, 'ss8'
    elif params_local['aa_ss'] == 'ss8' or params_local['aa_ss'] == 'ss3':
        print('--------- case 2: aa_ss = ss8 or ss3, PLEASE run FS_ss8_2_aa first --------------')

        train_data_ss8_df = pd.read_pickle(train_data_file)
        test_data_ss8_df = pd.read_pickle(test_data_file)
        SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'SPOT1DLM_aass3ss8.pkl')

        print('\ntrain_data_ss8_df = \n', train_data_ss8_df['proteins'])
        print('\ntest_data_ss8_df = \n', test_data_ss8_df['proteins'])
        print('\nSPOT1DLM_aass3ss8 = \n', SPOT1DLM_aass3ss8)

        SPOT1DLM_aass3ss8_dict = SPOT1DLM_aass3ss8.set_index('proteins').T.to_dict('list')


        def ss8_2_aa(df):
            k = 0
            for index, row in df.iterrows():
                # print(k)
                k += 1
                aa_replace = SPOT1DLM_aass3ss8_dict[row['proteins']][0]  # [0]:AA, [1]:SS3, [2]:SS8
                df.loc[index, 'sequences'] = aa_replace
            temp_df = df
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


        # pkl2fa
        print('pkl-2-fa starting ...')
        pkl2fa('data/train_data_aa.pkl', 'data/train_data.fa')
        pkl2fa('data/test_data_aa.pkl', 'data/test_data.fa')

    else:
        print('Holy crap! No specific cases...')

    print('\n################## And they all lived happily ever after! ##################\n')




if __name__ == '__main__':
    main()



