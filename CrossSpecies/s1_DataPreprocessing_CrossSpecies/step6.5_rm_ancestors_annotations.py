
import pandas as pd
from step0_DataPreprocessingSetting import *
from utils import Ontology, is_exp_code, is_cafa_target, FUNC_DICT


go_file = path_pub_data+'go.obo'
fpath_input = path_pub_data + 'swissprot_clean_pkl_NO_rm_annotations/swissprot_clean_ALL00_aa.pkl'
fpath_output = path_pub_data + 'swissprot_clean_ALL00_aa.pkl'

go = Ontology(go_file, with_rels=True)
df = pd.read_pickle(fpath_input)

print(df.info())


rm_ancestors_annotations = []
for i, row in df.iterrows():
    exp_annots = row['exp_annotations']
    rest_set = go.rm_ancestors(set(exp_annots))
    rest_list = list(rest_set)
    rm_ancestors_annotations.append(rest_list)
df['rm_ancestors_annotations'] = rm_ancestors_annotations

df.to_pickle(fpath_output)

print(print(df.info()))
print(df['rm_ancestors_annotations'])

