
import pandas as pd
import torch
import argparse
import numpy as np
# from tqdm import tqdm
# from dataset.data_functions import read_list, read_fasta_file  # original
# from data_functions import read_list, read_fasta_file
from tape import ProteinBertModel, TAPETokenizer
# from tape.models import ProteinBertModel
import esm
from esm.pretrained import load_model_and_alphabet_local



from step0_DataPreprocessingSetting import *




### FS ###
parser = argparse.ArgumentParser()
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')  # original: cpu

args = parser.parse_args()

### FS ###
model_path = '/home/fsong/work/py_proj/prot_data/esm/'  # esm1b_t33_650M_UR50S.pt
device = "cuda:0"  # cpu or "cuda:0" if you have a GPU

#
# tokenizer = TAPETokenizer()
# alphabet = tokenizer.vocab

model, alphabet = load_model_and_alphabet_local(model_path)

model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
batch_converter = alphabet.get_batch_converter()



# model, alphabet = torch.hub.load("facebookresearch/esm:main", "esm1b_t33_650M_UR50S")

# Load the model
# model = ProteinBertModel.from_pretrained(model_path).to('cuda:0')

# Load the tokenizer and get the alphabet


print('lol')
