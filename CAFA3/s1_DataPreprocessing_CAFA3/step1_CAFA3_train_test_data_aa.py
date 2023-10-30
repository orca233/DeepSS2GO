#!/usr/bin/env python



'''
要准备的原始数据就在这三个文件夹中：pub_data/data_cafa3/ + 
benchmark20171115/
CAFA3_targets/
CAFA3_training_data/


输入 INPUT：
1）train_sequences_file
uniprot_sprot_exp.fasta::::::::  66841 个文件
>A0A060X6Z0
MPISSSSSSSTKSMRRAASELERSD。。。。

2）train_annotations_file
uniprot_sprot_exp.txt：：：：：：：       
P48347	GO:0005524	F
P48347	GO:0009742	P
P48347	GO:0009737	P
Q9S9Z8	GO:0051117	F
...

3）test_sequences_file
targets_all.fasta：：：：：：： 130827 个文件？？？？
>T100900000001 1433B_MOUSE
MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTE。。。。

4）test_annotations_file
/benchmark20171115/groundtruth/leafonly_all.txt  ？？？？？
T100900005305	GO:0033234
T100900010085	GO:0006468
T100900010085	GO:0046777
T100900003996	GO:0001782


原始 deepgoplus 输出 OUTPUT：
out_terms_file  -- 4901 -- 'data_cafa3/terms.pkl'
train_data_file  -- 66841 -- 'data-cafa3/train_data.pkl'
test_data_file  -- 3328 -- 'data-cafa3/test_data.pkl'

更新 DeepSS2GO 输出 OUTPUT：
out_terms_file  -- 4901 -- 'data_cafa3/terms_all.pkl'  -- 没变 ---- 修改了，不在s1中生成terms了
train_data_file  -- 66841 -- 'data_cafa3/CAFA3_train_data_raw.pkl' -- 更新
test_data_file  -- 3328 -- 'data_cafa3/CAFA3_test_data_raw.pkl'  -- 更新

'''

import click as ck
import numpy as np
import pandas as pd
from collections import Counter
from utils import Ontology, FUNC_DICT, read_fasta
import logging
from step0_DataPreprocessingSetting import *

logging.basicConfig(level=logging.INFO)



@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')  # FS 添加
# @ck.option('--go-file', '-gf', default=path_base + 'pub_data/go.obo', help='Gene Ontology file in OBO Format')  # original: 'data_cafa3/go.obo'
### input
@ck.option('--train-sequences-file', '-trsf', default=path_pub_data + 'data_cafa3/CAFA3_training_data/uniprot_sprot_exp.fasta', help='CAFA training sequences fasta')
@ck.option('--train-annotations-file', '-traf', default=path_pub_data + 'data_cafa3/CAFA3_training_data/uniprot_sprot_exp.txt', help='CAFA training annotations fasta')
@ck.option('--test-sequences-file', '-tssf', default=path_pub_data + 'data_cafa3/CAFA3_targets/targets_all.fasta', help='CAFA training annotations fasta')
@ck.option('--test-annotations-file', '-tsaf', default=path_pub_data + 'data_cafa3/benchmark20171115/groundtruth/leafonly_all.txt', help='CAFA training annotations fasta')
## output:
@ck.option('--out-terms-file', '-otf', default=path_pub_data + 'data_cafa3/terms_all.pkl', help='Result file with a list of terms for prediction task')
@ck.option('--train-data-file', '-trdf', default=path_pub_data + 'data_cafa3/CAFA3_train_data_raw.pkl', help='Result file with a list of terms for prediction task')
@ck.option('--test-data-file', '-tsdf', default=path_pub_data + 'data_cafa3/CAFA3_test_data_raw.pkl', help='Result file with a list of terms for prediction task')
@ck.option('--min-count', '-mc', default=50, help='Minimum number of annotated proteins')

def main(go_file, train_sequences_file, train_annotations_file, test_sequences_file, test_annotations_file,
         out_terms_file, train_data_file, test_data_file, min_count):
    logging.info('Loading GO')
    go = Ontology(go_file, with_rels=True)
    
    logging.info('Loading training annotations')
    train_annots = {}
    with open(train_annotations_file, 'r') as f:
        for line in f:
            it = line.strip().split('\t')
            prot_id = it[0]
            if prot_id not in train_annots:
                train_annots[prot_id] = set()
            go_id = it[1]
            train_annots[prot_id].add(go_id)
                
    logging.info('Loading training sequences')
    info, seqs = read_fasta(train_sequences_file)  # info = A0A060X6Z0, seqs = MPISSSSSSSTKSMRRAAS...
    proteins = []
    sequences = []
    annotations = []

    for prot_info, sequence in zip(info, seqs):
        prot_id = prot_info.split()[0]
        if prot_id in train_annots:
            proteins.append(prot_id)
            sequences.append(sequence)
            annotations.append(train_annots[prot_id])
    
    prop_annotations = []
    cnt = Counter()
    for annots in annotations:
        # Propagate annotations
        annots_set = set()
        for go_id in annots:
            annots_set |= go.get_ancestors(go_id)  # original: go.get_anchestors(go_id)
        prop_annotations.append(annots_set)
        for go_id in annots_set:
            cnt[go_id] += 1
        
    df = pd.DataFrame({
        'proteins': proteins,
        'sequences': sequences,
        'prop_annotations': prop_annotations,  # original: 'annotations': prop_annotations,
    })

    logging.info(f'Train proteins: {len(df)}')
    logging.info(f'Saving training data to {train_data_file}')
    df.to_pickle(train_data_file)



    # 下面这段注释掉了，两个原因
    # 1）我们的terms不止>50，还要 train_data X test_data 的GO terms
    # 2）这里的到的 train_data & test_data只是 初始值，后面还要step2进行筛选。所以把terms_all.pkl放在 s2_TrainTest里计算

    # # Filter terms with annotations more than min_count
    # res = {}
    # for key, val in cnt.items():
    #     if val >= min_count:  # 对于'训练集'中大于50次的
    #         ont = key.split(':')[0]
    #         if ont not in res:
    #             res[ont] = []
    #         res[ont].append(key)
    # terms = []
    # for key, val in res.items():
    #     terms += val
    #
    # logging.info(f'Number of terms {len(terms)}')
    # logging.info(f'Saving terms to {out_terms_file}')
    #
    # df = pd.DataFrame({'terms': terms})
    # df.to_pickle(out_terms_file)



    logging.info('Loading testing annotations')
    test_annots = {}
    with open(test_annotations_file, 'r') as f:
        for line in f:
            it = line.strip().split('\t')
            prot_id = it[0]
            if prot_id not in test_annots:
                test_annots[prot_id] = set()
            go_id = it[1]
            test_annots[prot_id].add(go_id)
                
    logging.info('Loading testing sequences')
    info, seqs = read_fasta(test_sequences_file)
    proteins = []
    sequences = []
    annotations = []
    for prot_info, sequence in zip(info, seqs):
        prot_id = prot_info.split()[0]
        if prot_id in test_annots:  # 不是所有seq的protein名称都在leafonly_all的清单（包含对应的GO）中。
            proteins.append(prot_id)
            sequences.append(sequence)
            annotations.append(test_annots[prot_id])
    
    prop_annotations = []
    for annots in annotations:
        # Propagate annotations
        annots_set = set()
        for go_id in annots:
            annots_set |= go.get_ancestors(go_id)
        prop_annotations.append(annots_set)
        
    df = pd.DataFrame({
        'proteins': proteins,
        'sequences': sequences,
        'prop_annotations': prop_annotations,  # original: 'annotations': prop_annotations
    })
    logging.info(f'Test proteins {len(df)}')
    logging.info(f'Saving testing data to {test_data_file}')
    df.to_pickle(test_data_file)

                


if __name__ == '__main__':
    main()
