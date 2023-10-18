
import pandas as pd
from step0_DataPreprocessingSetting import *


def pkl2fa(input_file, output_file):  # DeepGOPlus 原文件为 diamond_data.py
    # Load interpro data
    df = pd.read_pickle(input_file)
    print('pkl length', len(df))
    with open(output_file, 'w') as f:
        for row in df.itertuples():
            f.write('>' + row.proteins + '\n')
            f.write(row.sequences + '\n')


input_file = path_base + 'pub_data/data_new/new_clean_aa.pkl'
output_file = path_base + 'pub_data/data_new/new_clean_aa.fa'
pkl2fa(input_file, output_file)


input_file = path_base + 'pub_data/data_new/new_clean_ss8.pkl'
output_file = path_base + 'pub_data/data_new/new_clean_ss8.fa'
pkl2fa(input_file, output_file)


