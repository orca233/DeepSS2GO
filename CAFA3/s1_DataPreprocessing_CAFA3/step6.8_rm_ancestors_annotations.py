# README: 把swissprot_clean_ALL00_aa.pkl 添加一列，添加列为去掉所有 “祖先GO” rm_ancestors_annotations的。保留最孙子辈
# 输入文件为原先的在 swissprot_clean_pkl_NO_rm_annotations/ 文件夹下的 swissprot_clean_ALL00_aa.pkl
# 输出直接到 put_data/ 中

import pandas as pd
from step0_DataPreprocessingSetting import *
from utils import Ontology, is_exp_code, is_cafa_target, FUNC_DICT


go_file = path_pub_data+'go.obo'
fpath_input = path_pub_data + 'swissprot_clean_pkl_NO_rm_annotations/swissprot_clean_ALL00_aa.pkl'
fpath_output = path_pub_data + 'swissprot_clean_ALL00_aa.pkl'  # 添加新列 df['rm_ancestors_annotations']

go = Ontology(go_file, with_rels=True)
df = pd.read_pickle(fpath_input)

print(df.info())

######  FS add，除去所有父类，只保留最孙子的节点，
# 可以直接从exp_annotation里删就行，因为prop_annot是exp基础上并入的所有祖先 -- 类似于互质 ###
rm_ancestors_annotations = []  # prop: 父类的GO，往上找
for i, row in df.iterrows():
    # rm_ancestors annotations
    exp_annots = row['exp_annotations']
    rest_set = go.rm_ancestors(set(exp_annots))
    rest_list = list(rest_set)
    rm_ancestors_annotations.append(rest_list)
df['rm_ancestors_annotations'] = rm_ancestors_annotations
######  FS add, -- done ###

df.to_pickle(fpath_output)

#
print('111111')
print(print(df.info()))
print(df['rm_ancestors_annotations'])

