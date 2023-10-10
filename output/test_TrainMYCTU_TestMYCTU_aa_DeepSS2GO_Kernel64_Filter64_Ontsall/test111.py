import pandas as pd
path_pub_data = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/pub_data/'

SPOT1DLM_aass3ss8 = pd.read_pickle(path_pub_data + 'SPOT1DLM_aass3ss8.pkl')

print(SPOT1DLM_aass3ss8.info())

print('2222222222222222')
fpath = '~/work/py_proj/prot_algo/DeepSS2GO_AB/output/test_TrainMYCTU_TestMYCTU_ss8_DeepSS2GO_Kernel128_Filter128_Ontsall/data/test_data_aa.pkl'
df = pd.read_pickle(fpath)
print(df.info())
print(df['sequences'])
