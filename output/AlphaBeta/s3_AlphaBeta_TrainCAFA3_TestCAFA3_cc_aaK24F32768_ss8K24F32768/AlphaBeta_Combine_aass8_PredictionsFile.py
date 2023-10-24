# README 这个文件可以把predictions_aa.pkl的preds这一列重命名为preds_aa
# 把predictions_ss8.pkl的preds这一列重命名为preds_ss8
# 把 preds_aa 添到 predictions_ss8.pkl 的后面

import pandas as pd


# fpath = 'data/predictions_ss8.pkl'

df_aa = pd.read_pickle('data/predictions_aa.pkl')
df_ss8 = pd.read_pickle('data/predictions_ss8.pkl')

# 将"preds"列重命名为"preds_aa"
df_aa.rename(columns={'preds': 'preds_aa'}, inplace=True)
df_ss8.rename(columns={'preds': 'preds_ss8'}, inplace=True)

# 将重命名后的"preds_ss8"列添加到df_aa的最后一列
df_aa['preds_ss8'] = df_ss8['preds_ss8']

# 输出结果
print(df_aa)
print(df_aa.info())

df_aa.to_pickle('data/predictions_aa_ss8.pkl')


print('happy ending')
