
import pandas as pd
import torch
import argparse
import numpy as np
import esm
from step0_DataPreprocessingSetting import *

save_path_npy = path_base + 'redundancy/SPOT1DLM_inputs_new/'
os.system('mkdir -p %s' % save_path_npy)
parser = argparse.ArgumentParser()
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')
args = parser.parse_args()
model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()
model = model.to(args.device)

fpath = path_base + 'pub_data/data_new/new_clean_aa.pkl'
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



