
import time
import os.path


import pandas as pd
import torch
import argparse
import numpy as np
# from tqdm import tqdm
# from dataset.data_functions import read_list, read_fasta_file  # original
# from data_functions import read_list, read_fasta_file
# from tape import ProteinBertModel, TAPETokenizer

# import esm
# from esm.pretrained import load_model_and_alphabet_local

from step0_DataPreprocessingSetting import *




### FS ###
parser = argparse.ArgumentParser()
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')  # original: cpu

args = parser.parse_args()

###########################################
##### 下载 esm1b_t33_650M_UR50S 三种方法 ####
##########################################

### 方法一：original: load model  适合在SUSTech上运行 ###
model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")
# # Using cache found in /home/songfu/.cache/torch/hub/facebookresearch_esm_main  在线找数据？啥也没有

### 方法二： FS try 2 用esm方法下载： 把hub下载好的加载到对应文件夹，适合在lab——linux_3090上运行 #####
# model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
# /home/fsong/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt & esm1b_t33_650M_UR50S.pt
### FS try 2 done ###


# ### 方法三 FS: load model  不好使。。。 ###
# model_path = '/home/fsong/work/py_proj/prot_data/esm/esm1b_t33_650M_UR50S.pt'
# device = "cuda:0"  # cpu or "cuda:0" if you have a GPU
# model = torch.load(model_path)  #.to(device)
# tokenizer = TAPETokenizer()
# alphabet = tokenizer.vocab
# ### FS change done ###






batch_converter = alphabet.get_batch_converter()
model = model.to(args.device)

# original:'swissprot_clean_ALL00_aa.pkl ------ 已根据CAFA3更新：------
fpath = path_pub_data + 'data_cafa3/CAFA3_test_data_clean_aa.pkl'  # new1: CAFA3_clean_train_data_aa.pkl   new2: CAFA3_clean_test_data_aa.pkl  修改！！！！！！！！！！！！！！
prot_df = pd.read_pickle(fpath)  # 包含prot_name & seq

esm_count = 0  # 统计生成了多少个esm文件，  i: 第几轮，esm_count:本轮进行了多少个
for index, row in prot_df.iterrows():
    prot_name = row['proteins']
    # save_path = save_path_npy + prot_name + "_esm.npy"  # 文件太大，存在redundancy里  # 这个是通用的swissprot，好像CAFA3的train&test命名应该不同把。。。我猜的。。。
    save_path = path_base + 'redundancy/SPOT1DLM_inputs_CAFA3_test_data/' + prot_name + "_esm.npy"  # 文件太大，存在redundancy里  # 这个针对CAFA3，因为train/test已经分化

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





#### original: #####

# for prot_path in tqdm(prot_list):
#
#     prot_name = prot_path.split('/')[-1].split('.')[0]
#     save_path = "inputs/" + prot_name + "_esm.npy"
#     seq = read_fasta_file(prot_path)
#     data = [(prot_name, seq)]
#
#     batch_labels, batch_strs, batch_tokens = batch_converter(data)
#     batch_tokens = batch_tokens.to(args.device)
#     with torch.no_grad():
#         results = model(batch_tokens, repr_layers=[33], return_contacts=True)
#     token_representations = results["representations"][33]
#
#     sequence_representations = []
#     for i, (prot_n, seq) in enumerate(data):
#         if args.device == "cpu":
#             save_arr = token_representations[i, 1: len(seq) + 1].numpy()
#         else:
#             save_arr = token_representations[i, 1: len(seq) + 1].cpu().numpy()
#         np.save(save_path, save_arr)
# print(" ESM-1b embeddings generation completed ... ")