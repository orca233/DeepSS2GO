import time
start_time = time.time()
import re
import torch
import argparse
import numpy as np
from transformers import T5EncoderModel, T5Tokenizer
from step0_DataPreprocessingSetting import *
import pandas as pd

# # step3已经建立了 redundancy/SPOT1DLM_inputs/
# # 用来储存esm/prottrans.npy文件/scem/work/songfu/prot_algo/DeepSS2GO/redundancy/sSPOT1DLM_inputs
# os.system('mkdir -p %sSPOT1DLM_inputs' % path_redundancy)

print('\n################## a long, long time ago ... ##################\n')
print('# starting step4_SPOT1DLM_generate_prottrans #')

parser = argparse.ArgumentParser()
# parser.add_argument('--file_list', default='', type=str, help='file list path ')
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')  # original: cpu
args = parser.parse_args()

### original:
# tokenizer = T5Tokenizer.from_pretrained("Rostlab/prot_t5_xl_uniref50", do_lower_case=False)
# model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_uniref50")

### FS 提前下载好 '/scem/work/songfu/prot_data/Prot_T5_XL_UniRef50'
tokenizer = T5Tokenizer.from_pretrained(path_Prot_T5_XL_UniRef50, do_lower_case=False)  #
model = T5EncoderModel.from_pretrained(path_Prot_T5_XL_UniRef50)

device = torch.device(args.device)
model = model.to(args.device)
model = model.eval()


### FS ###
prot_df = pd.read_pickle(path_pub_data + 'swissprot_clean_ALL00_aa.pkl')

k = 0  # 统计个数
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

print('\n################## And they all lived happily ever after! ##################\n')
