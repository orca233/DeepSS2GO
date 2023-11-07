##!/usr/bin/env python

import click as ck
import numpy as np
import pandas as pd
import gzip

from collections import Counter
# from aminoacids import MAXLEN, to_ngrams
import logging

logging.basicConfig(level=logging.INFO)

# DeepGOPlus 原文件为 diamond_data.py

# README: 这是一个pkl2pkl程序，把train_data.pkl中的10-12列，只提取prop_annotations这一列

@ck.command()
@ck.option('--data-file', '-df', default='data/train_data_backup.pkl', help='Pandas dataframe with protein sequences')
@ck.option('--out-file', '-o', default='data/train_data.pkl', help='Fasta file')

def main(data_file, out_file):
    # Load interpro data
    df = pd.read_pickle(data_file)
    print('pkl length', len(df))

    df = df.drop(columns=['index', 'accessions', 'sequences', 'annotations', 'interpros', 'orgs', 'exp_annotations', 'cafa_target'])

    print(df.info())
    print(df)
    df.to_pickle(out_file)

if __name__ == '__main__':
    main()
