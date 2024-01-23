import pandas as pd
import torch
import argparse
import numpy as np
from step0_DataPreprocessingSetting import *

print('\n################## a long, long time ago ... ##################\n')
print('# starting step3_SPOT1DLM_generage_esm #')


parser = argparse.ArgumentParser()
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')

args = parser.parse_args()



model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")


batch_converter = alphabet.get_batch_converter()
model = model.to(args.device)

fpath = path_base + 'pub_data/swissprot_clean_ALL00_aa.pkl'
prot_df = pd.read_pickle(fpath)

esm_count = 0
for index, row in prot_df.iterrows():
    prot_name = row['proteins']
    save_path = save_path_npy + prot_name + "_esm.npy"
    seq = row['sequences']
    data = [(prot_name, seq)]

    print(fpath + '-' + str(esm_count))
    esm_count += 1
    print(prot_name)

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_tokens = batch_tokens.to(args.device)
    with torch.no_grad():
        results = model(batch_tokens, repr_layers=[33], return_contacts=True)
    token_representations = results["representations"][33]

    sequence_representations = []
    for i, (prot_n, seq) in enumerate(data):
        if args.device == "cpu":
            save_arr = token_representations[i, 1: len(seq) + 1].numpy()
        else:   # GPU
            save_arr = token_representations[i, 1: len(seq) + 1].cpu().numpy()
        np.save(save_path, save_arr)

print(" ESM-1b embeddings generation completed ... ")

print('\n################## And they all lived happily ever after! ##################\n')

