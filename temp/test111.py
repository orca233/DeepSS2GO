import pandas as pd
import numpy as np

fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/'

df = pd.read_pickle(fpath + 'swissprot_clean_ALL00_aa.pkl')

print(df.info())
