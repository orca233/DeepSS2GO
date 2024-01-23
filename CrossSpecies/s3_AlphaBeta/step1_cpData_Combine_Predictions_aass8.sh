##!/bin/bash

mkdir -p data
path_current="$(pwd)/"

path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO/"


path_aa="${path_base}output/Best4PredictNew/test_TrainALL00_TestALL00_aa_DeepSS2GO_Kernel16_Filter65536_Ontsall/"
path_ss8="${path_base}output/Best4PredictNew/test_TrainALL00_TestALL00_ss8_DeepSS2GO_Kernel48_Filter8192_Ontsall/"



cp "${path_aa}data/go.obo" "${path_current}data/go.obo"

# aa
cp "${path_aa}data/train_data.pkl" "${path_current}data/train_data_aa.pkl"
cp "${path_aa}data/train_data.fa" "${path_current}data/train_data_aa.fa"

cp "${path_aa}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_aa.pkl"
cp "${path_aa}data/predictions.pkl" "${path_current}data/predictions_aa.pkl"
cp "${path_aa}step0_TrainTestSetting_local.py" "${path_current}step0_TrainTestSetting_local_aa.py"
cp "${path_aa}data/test_data.pkl" "${path_current}data/test_data_aa.pkl"
cp "${path_aa}data/test_data.fa" "${path_current}data/test_data_aa.fa"
cp "${path_aa}data/model_checkpoint.pth" "${path_current}data/model_checkpoint_aa.pth"

# ss8
cp "${path_ss8}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_ss8.pkl"
cp "${path_ss8}data/predictions.pkl" "${path_current}data/predictions_ss8.pkl"
cp "${path_ss8}step0_TrainTestSetting_local.py" "${path_current}step0_TrainTestSetting_local_ss8.py"
cp "${path_ss8}data/test_data.pkl" "${path_current}data/test_data_ss8.pkl"
cp "${path_ss8}data/test_data.fa" "${path_current}data/test_data_ss8.fa"
cp "${path_ss8}data/model_checkpoint.pth" "${path_current}data/model_checkpoint_ss8.pth"


python AlphaBeta_Combine_aass8_PredictionsFile.py


