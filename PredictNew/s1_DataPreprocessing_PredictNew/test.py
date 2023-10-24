import pandas as pd

fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/data_new/new_clean_aa.pkl'
df = pd.read_pickle(fpath)

print(df.info())

