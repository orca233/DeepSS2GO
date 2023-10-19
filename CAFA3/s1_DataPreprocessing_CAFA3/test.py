import pandas as pd

fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/pub_data/data_cafa3/'

df = pd.read_pickle(fpath + 'test_data.pkl')
print(df.info())

print(df.info)
