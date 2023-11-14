
# README
# 把esm & prottrans的结果整合，*esm.npy & *pt.npy 各有71000个文件
# 联合SPOT1DLM自己算的checkpoint得到新的输出结果，结果包含ss3, ss8和溶解度等信息


import torch
import argparse
from torch.utils.data import DataLoader
# from config import PATH, LIST, VAL_LIST, TEST_LIST, TEST2_LIST, TEST3_LIST, TEST4_LIST, IGNORE_LABEL, DEVICE
import pandas as pd
from dataset_inference import text_collate_fn, FS_Proteins_Dataset  # FS_Proteins_Dataset会调用inputs/*esm.npy *pt.npy

from models.bilstm import Network
from models.ms_resnet import Network as Network2
from models.ms_res_lstm import Network as Network3
from models.bilstm_reg import Network as Network4
from models.ms_resnet_reg import Network as Network5
from models.ms_res_lstm_reg import Network as Network6

from main import main_reg, main_class, write_csv

from step0_DataPreprocessingSetting import *

save_path_csv = path_redundancy + 'SPOT1DLM_results_new/'  # for step5 inference生成的 *csv文件
os.system('mkdir -p %s' % save_path_csv)

parser = argparse.ArgumentParser()
parser.add_argument('--file_list', default='', type=str, help='file list path ')
parser.add_argument('--device', default='cuda:0', type=str, help=' define the device you want the ')  # original: cpu
parser.add_argument('--save_path', default=save_path_csv, type=str, help=' define the device you want the ')
args = parser.parse_args()


### FS ###
prot_df = pd.read_pickle(path_base + 'pub_data/data_new/new_clean_aa.pkl')

FS_test_set = FS_Proteins_Dataset(prot_df)  ## FS: 俺也不知道好使不。。。 看起来还OJBK
print("test_dataset Loaded with ", len(FS_test_set), "proteins")
# this implementation has only been tested for batch size 1 only.
test_loader = DataLoader(FS_test_set, batch_size=1, collate_fn=text_collate_fn, num_workers=16)


# original
# test_set = Proteins_Dataset(args.file_list)  ## spot-1d test set
# print("test_dataset Loaded with ", len(test_set), "proteins")
# # this implementation has only been tested for batch size 1 only.
# test_loader = DataLoader(test_set, batch_size=1, collate_fn=text_collate_fn, num_workers=16)


torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

model1 = Network()
model2 = Network2()
model3 = Network3()
model4 = Network4()
model5 = Network5()
model6 = Network6()

model1.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model1.pt", map_location=torch.device('cpu')))  # original: device('cpu') 换成cuda:0不好使
model2.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model2.pt", map_location=torch.device('cpu')))
model3.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model3.pt", map_location=torch.device('cpu')))
model4.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model4.pt", map_location=torch.device('cpu')))
model5.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model5.pt", map_location=torch.device('cpu')))
model6.load_state_dict(torch.load(path_SPOT1DLM_checkpoints + "model6.pt", map_location=torch.device('cpu')))

model1 = model1.to(args.device)
model2 = model2.to(args.device)
model3 = model3.to(args.device)
model4 = model4.to(args.device)
model5 = model5.to(args.device)
model6 = model6.to(args.device)


class_out = main_class(test_loader, model1, model2, model3, args.device)
names, seq, ss3_pred_list, ss8_pred_list, ss3_prob_list, ss8_prob_list = class_out
reg_out = main_reg(test_loader, model4, model5, model6, args.device)
psi_list, phi_list, theta_list, tau_list, hseu_list, hsed_list, cn_list, asa_list = reg_out
print(len(ss3_pred_list), len(psi_list))
write_csv(class_out, reg_out, args.save_path)  # save

