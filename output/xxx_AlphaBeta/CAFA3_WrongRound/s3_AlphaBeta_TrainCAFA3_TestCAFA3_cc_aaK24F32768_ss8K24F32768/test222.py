#!/usr/bin/env python


import numpy as np
import pandas as pd
import click as ck
# import logging
# from utils import FUNC_DICT, Ontology, NAMESPACES
# from step0_TrainTestSetting_global import *


@ck.command()
@ck.option('--terms-file', '-tf', default='./data/terms_gominre_trxte_aa.pkl', help='set of annotations')  # original 'data/terms_all.pkl'



def main(terms_file):  # FS加的go_file

    terms_df_aa = pd.read_pickle(terms_file)
    print(terms_df_aa)








#####################################################################################

if __name__ == '__main__':
    alphas = {'mf': 0, 'bp': 0, 'cc': 0}
    main()
