import pandas as pd
import numpy as np


fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/swissprot.pkl'

df = pd.read_pickle(fpath)

print(type(df['proteins']))


for i in df['proteins']:
    if 'NEIMB' in i:
        print(i)


