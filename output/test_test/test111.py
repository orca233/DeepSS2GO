import pandas as pd
# path_pub_data = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/pub_data/'
#
# SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'SPOT1DLM_aass3ss8.pkl')
#
# print(SPOT1DLM_aass3ss8.info())

print('2222222222222222')
fpath = '/data0/fsong/py_proj/prot_algo/DeepSS2GO_AB/output/test_test/data/new_aa.pkl'
df = pd.read_pickle(fpath)
print(df.info())
# print(df['sequences'])

print(df['proteins'])



