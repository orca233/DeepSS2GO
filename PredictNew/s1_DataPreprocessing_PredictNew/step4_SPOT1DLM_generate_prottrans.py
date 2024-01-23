import time
start_time = time.time()
import re
import torch
import argparse
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from step0_DataPreprocessingSetting import *
import pandas as pd


save_path_npy = path_base + 'redundancy/SPOT1DLM_inputs_new/'
os.system('mkdir -p %s' % save_path_npy)

parser = argparse.ArgumentParser()
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')
args = parser.parse_args()

tokenizer = T5Tokenizer.from_pretrained(path_Prot_T5_XL_UniRef50, do_lower_case=False)  #
model = T5EncoderModel.from_pretrained(path_Prot_T5_XL_UniRef50)

device = torch.device(args.device)
model = model.to(args.device)
model = model.eval()



prot_df = pd.read_pickle(path_base + 'pub_data/data_new/new_clean_aa.pkl')

k = 0
for index, row in prot_df.iterrows():
    print(k)
    k += 1

    seq = row['sequences']
    prot_name = row['proteins']

    seq_temp = seq.replace('', " ")
    sequences_Example = [seq_temp]
    sequences_Example = [re.sub(r"[UZOB]", "X", sequence) for sequence in sequences_Example]
    ids = tokenizer.batch_encode_plus(sequences_Example, add_special_tokens=True, padding=True)

    input_ids = torch.tensor(ids['input_ids']).to(args.device)
    attention_mask = torch.tensor(ids['attention_mask']).to(args.device)
    with torch.no_grad():
        embedding = model(input_ids=input_ids, attention_mask=attention_mask)

    if args.device == "cpu":
        embedding = embedding.last_hidden_state.numpy()
    else:
        embedding = embedding.last_hidden_state.cpu().numpy()

    features = []
    for seq_num in range(len(embedding)):
        seq_len = (attention_mask[seq_num] == 1).sum()
        seq_emd = embedding[seq_num][:seq_len - 1]
        features.append(seq_emd)

    np.save(save_path_npy + prot_name + "_pt.npy", features[0])


print(" ProtTrans embeddings generation completed ... ")
