##!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
from collections import Counter
import logging
import os
from utils import Ontology, NAMESPACES
from step0_TrainTestSetting_local import *
from step0_TrainTestSetting_global import path_base

logging.basicConfig(level=logging.INFO)


os.system('mkdir -p data')
os.system('mkdir -p results')

@ck.command()
@ck.option('--go-file', '-gf', default='data/go.obo', help='Gene Ontology file in OBO Format')
@ck.option('--train-data-file', '-trdf', default='data/train_data.pkl', help='train data file')
@ck.option('--test-data-file', '-tedf', default='data/test_data.pkl', help='XX')  # output
@ck.option('--terms-gominrepeat-file', '-tgf', default='data/terms_gominre.pkl', help='XX')
@ck.option('--terms-gominrepeat-trainxtest-file', '-tgtf', default='data/terms_gominre_trxte.pkl', help='XX')
@ck.option('--terms-bp-file', '-tbf', default='data/terms_bp.pkl', help='XX')
@ck.option('--terms-cc-file', '-tcf', default='data/terms_cc.pkl', help='XX')
@ck.option('--terms-mf-file', '-tmf', default='data/terms_mf.pkl', help='XX')
@ck.option('--train-data', '-trd', default=params_local['train_data'], help='HUMAN/ECOLI/...')
@ck.option('--test-data', '-ted', default=params_local['test_data'], help='HUMAN/ECOLI/...')
@ck.option('--aa-ss', '-aass', default=params_local['aa_ss'], help='aa/ss8/ss3')
@ck.option('--gominrepeat', '-gmr', default=params_local['GOMinRepeat'], help='GO min repeat times')
@ck.option('--path-base', '-pb', default=params_local['path_base'], help='..')


def main(go_file, train_data_file, test_data_file, terms_gominrepeat_file, terms_gominrepeat_trainxtest_file,
         terms_bp_file, terms_cc_file, terms_mf_file, train_data, test_data, aa_ss, gominrepeat, path_base):
    print('\n################## a long, long time ago ... ##################\n')
    print('# starting step1_SplitTrainTest_Terms #')
    print('# params_local:')
    for key, value in params_local.items():
        print(key, ' = ', value)


    ######### pre-setting ###################

    path_pub_data = path_base + 'pub_data/'
    path_redundancy = path_base + 'redundancy/'



    ###### step1.1  ######
    print('\n-------- Step 1.1: Processing train_data_file & test_data_file ------------')
    if train_data == test_data and train_data != 'CAFA3':
        print('\n----- case 0 :train_data == test_data != CAFA3 --------')
        swissprot_data_file = path_pub_data + 'swissprot_clean_%s_%s.pkl' % (train_data, aa_ss)
        all_data_df = pd.read_pickle(swissprot_data_file)

        n = len(all_data_df)
        index = np.arange(n)
        train_n = int(n * params_local['TrainTestRatio'])
        np.random.seed(seed=0)  # ？
        np.random.shuffle(index)
        train_df = all_data_df.iloc[index[:train_n]]
        test_df = all_data_df.iloc[index[train_n:]]

        # save pkl
        train_df.to_pickle(train_data_file)
        test_df.to_pickle(test_data_file)

        print('train data = test data = %s' % train_data)
        print("Number of all proteins", len(all_data_df))
        print('Number of train proteins', len(train_df))
        print('Number of test proteins', len(test_df))


    elif train_data != test_data and train_data != 'CAFA3' and test_data != 'CAFA3':
        print('\n----- case 1: train_data !!!=== test_data ------')  # logging.info
        os.system('cp %sswissprot_clean_%s_%s.pkl data/train_data.pkl' % (path_pub_data, train_data, aa_ss))
        os.system('cp %sswissprot_clean_%s_%s.pkl data/test_data.pkl' % (path_pub_data, test_data, aa_ss))

        # read pkl
        train_df = pd.read_pickle('data/train_data.pkl')
        test_df = pd.read_pickle('data/test_data.pkl')
        all_data_df = pd.concat([train_df, test_df]).reset_index(drop=True)

        print('train data = %s' % train_data)
        print('test data = %s' % test_data)
        print("Number of all proteins (containing 2 species)", len(all_data_df))
        print('Number of train proteins', len(train_df))
        print('Number of test proteins', len(test_df))


    elif train_data == test_data and train_data == 'CAFA3':
        print('\n------ case 2: train_data == test_data == CAFA3 -------')
        os.system('cp %sdata_cafa3/CAFA3_train_data_clean_%s.pkl data/train_data.pkl' % (path_pub_data, aa_ss))
        os.system('cp %sdata_cafa3/CAFA3_test_data_clean_%s.pkl data/test_data.pkl' % (path_pub_data, aa_ss))

        # read pkl
        train_df = pd.read_pickle('data/train_data.pkl')
        test_df = pd.read_pickle('data/test_data.pkl')
        all_data_df = pd.concat([train_df, test_df]).reset_index(drop=True)

        print('train data = %s' % train_data)
        print('test data = %s' % test_data)
        print("Number of all proteins (containing 2 species)", len(all_data_df))
        print('Number of train proteins', len(train_df))
        print('Number of test proteins', len(test_df))



    ########## step1.2  ###############

    print('\n-----Step 1.2: --------- Processing terms ---------------------')  # logging.info
    print('GO Min Repeat=', gominrepeat)
    cnt = Counter()
    for i, row in all_data_df.iterrows():
        for term in row['prop_annotations']:
            cnt[term] += 1

    print('Number of terms (raw total = train + test) = ', len(cnt))

    # Filter terms with annotations more than gominre
    res = {}
    for key, val in cnt.items():
        if val >= gominrepeat:
            ont = key.split(':')[0]
            if ont not in res:
                res[ont] = []
            res[ont].append(key)

    terms = []
    for key, val in res.items():
        terms += val

    terms_GOMinRepeat = terms
    print(f'Number of terms (> gominre) = : {len(terms_GOMinRepeat)}')


    train_prop_annotation = set()
    test_prop_annotation = set()
    for index, row in train_df.iterrows():
        prop_annotations = train_df.loc[index, 'prop_annotations']
        train_prop_annotation = train_prop_annotation | set(prop_annotations)
    for index, row in test_df.iterrows():
        prop_annotations = test_df.loc[index, 'prop_annotations']  #
        test_prop_annotation = test_prop_annotation | set(prop_annotations)

    train_x_test_prop_annotations = train_prop_annotation & test_prop_annotation
    terms_trainXtest = train_x_test_prop_annotations
    print('Number of terms (trainXtest) =', len(terms_trainXtest))

    terms_GOMinRe_trXte = terms_trainXtest & set(terms_GOMinRepeat)
    print('Number of terms ( >gominre & train_x_test_terms ) = ', len(terms_GOMinRe_trXte))


    # Save
    terms_GOMinRepeat = set(terms_GOMinRepeat)
    terms_df_GOMinRepeat = pd.DataFrame({'terms': list(terms_GOMinRepeat)})
    terms_df_GOMinRepeat.to_pickle(terms_gominrepeat_file)

    terms_df_gominre_trxte = pd.DataFrame({'terms': list(terms_GOMinRe_trXte)})
    terms_df_gominre_trxte.to_pickle(terms_gominrepeat_trainxtest_file)


    go_rels = Ontology(go_file, with_rels=True)

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['bp'])
    terms_bp_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_bp_df.to_pickle(terms_bp_file)
    print('terms_bp len=', len(terms_bp_df))  # 1517

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['cc'])
    terms_cc_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_cc_df.to_pickle(terms_cc_file)
    print('terms_cc len=', len(terms_cc_df))  # 269

    go_set_ont = go_rels.get_namespace_terms(NAMESPACES['mf'])
    terms_mf_df = terms_df_gominre_trxte[terms_df_gominre_trxte['terms'].isin(go_set_ont)]
    terms_mf_df.to_pickle(terms_mf_file)
    print('terms_mf len=', len(terms_mf_df))  # 302



    print('\n################## And they all lived happily ever after! ##################\n')

if __name__ == '__main__':
    main()






