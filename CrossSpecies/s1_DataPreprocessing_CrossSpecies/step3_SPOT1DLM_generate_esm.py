import pandas as pd
import torch
import argparse
import numpy as np
from step0_DataPreprocessingSetting import *

print('\n################## a long, long time ago ... ##################\n')
print('# starting step3_SPOT1DLM_generage_esm #')

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

# ### 方法三 FS: load model  不好使。。。 ###
# model_path = '/home/fsong/work/py_proj/prot_data/esm/esm1b_t33_650M_UR50S.pt'
# device = "cuda:0"  # cpu or "cuda:0" if you have a GPU
# model = torch.load(model_path)  #.to(device)
# tokenizer = TAPETokenizer()
# alphabet = tokenizer.vocab
# ### FS change done ###

batch_converter = alphabet.get_batch_converter()
model = model.to(args.device)

fpath = path_base + 'pub_data/swissprot_clean_ALL00_aa.pkl'  # 对于所有 swissprot_cleanup 一起算71000个
prot_df = pd.read_pickle(fpath)  # 包含prot_name & seq

esm_count = 0  # 统计生成了多少个esm文件，  i: 第几轮，esm_count:本轮进行了多少个
for index, row in prot_df.iterrows():
    prot_name = row['proteins']
    save_path = save_path_npy + prot_name + "_esm.npy"  # 文件太大，存在redundancy里
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

