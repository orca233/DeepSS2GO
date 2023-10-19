import pandas as pd
from step0_DataPreprocessingSetting import *



# output = []
# for keep_species in ['HUMAN', 'ARATH', 'MOUSE', 'YEAST', 'SCHPO', 'RAT', 'ECOLI', 'DROME', 'CAEEL', 'MYCTU']:
#     fpath = '%sswissprot_clean_%s_aa.pkl' % (path_pub_data, keep_species)
#     df = pd.read_pickle(fpath)
#     species_num = '%s: %d' % (keep_species, len(df))
#     output.append(species_num)
#
# print(output)


fpath = '/s0/fsong/py_proj/prot_algo/DeepSS2GO_Pytorch/pub_data/data_cafa3/CAFA3_train_data_raw.pkl'

df = pd.read_pickle(fpath)

print(df)
print(df.info())



