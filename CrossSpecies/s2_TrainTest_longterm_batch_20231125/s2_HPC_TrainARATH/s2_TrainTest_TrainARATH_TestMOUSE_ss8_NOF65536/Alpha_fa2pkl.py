##!/usr/bin/env python

import click as ck
import pandas as pd


@ck.command()
@ck.option('--input-file', '-inf', default='data/x.fa')
@ck.option('--output-file', '-ouf', default='data/x.pkl')

def main(input_file, output_file):
    sequences = read_fa_file(input_file)
    # print('111111111111111111111')
    # print(sequences)
    df = pd.DataFrame(sequences, columns=['proteins', 'sequences', 'prop_annotations'])  # original 'index', 'seq', 'prop_anno'
    df.to_pickle(output_file)
    print(f"Processed {len(sequences)} sequences and saved to {output_file}")
    print('233333333333333333333')



def read_fa_file(file_path):
    ind_sequences_prop = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        index = None
        sequence = ""
        prop_anno = ''  # 需要创建一个虚假的傀儡的prop_anno，否则step5_predict会报错
        for line in lines:
            if line.startswith('>'):
                if index is not None:
                    ind_sequences_prop.append((index, sequence, prop_anno))
                index = str(line[1:].strip())  # original: index = int(line[1:].strip())
                sequence = ""
            else:
                sequence += line.strip()
        if index is not None and sequence:
            ind_sequences_prop.append((index, sequence, prop_anno))

    return ind_sequences_prop




if __name__ == '__main__':
    main()


