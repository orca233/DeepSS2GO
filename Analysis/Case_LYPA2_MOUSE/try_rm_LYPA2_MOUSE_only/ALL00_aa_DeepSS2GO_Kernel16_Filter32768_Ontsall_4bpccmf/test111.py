import pandas as pd

fpath = './data/train_data.pkl'

df = pd.read_pickle(fpath)

print(df['proteins'])




