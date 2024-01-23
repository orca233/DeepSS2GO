import pandas as pd

df_aa = pd.read_pickle('data/predictions_aa.pkl')
df_ss8 = pd.read_pickle('data/predictions_ss8.pkl')

df_aa.rename(columns={'preds': 'preds_aa'}, inplace=True)
df_ss8.rename(columns={'preds': 'preds_ss8'}, inplace=True)

df_aa['preds_ss8'] = df_ss8['preds_ss8']

print(df_aa)
print(df_aa.info())

df_aa.to_pickle('data/predictions_aa_ss8.pkl')


print('happy ending')
