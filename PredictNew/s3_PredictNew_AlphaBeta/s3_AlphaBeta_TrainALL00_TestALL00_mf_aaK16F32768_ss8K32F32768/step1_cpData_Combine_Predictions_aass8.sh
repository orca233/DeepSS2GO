##!/bin/bash


# 当前文件夹要重新命名，i.e. step1_FindAlphaBeta_TrainMYCTU_TestMYCTU_aaK64F64_ss8K128F128


# 指定文件夹, aa & ss8
mkdir -p data
path_current="$(pwd)/"  # path_current后面没有/

# 绝对路径
path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO/"
# 相对路径
#path_base=$(cd ../../ && pwd)
#path_base="$(path_base)/"

path_aa="${path_base}output/test_TrainALL00/test_TrainALL00_TestALL00_aa/DeepSS2GO_Kernel16_Filter32768_Ontsall/"
path_ss8="${path_base}output/test_TrainALL00/test_TrainALL00_TestALL00_ss8/DeepSS2GO_Kernel32_Filter32768_Ontsall/"


###################### predictions.pkl & model_checkpoint.pth 的aa/ss8 最好还是用 link -s .... 要不太大了 #############

#### cp ####
# go.obo
cp "${path_aa}data/go.obo" "${path_current}data/go.obo"

# aa
#cp "${path_aa}data/diamond_aa.res" "${path_current}data/diamond_aa.res"
cp "${path_aa}data/train_data.pkl" "${path_current}data/train_data_aa.pkl"
cp "${path_aa}data/train_data.fa" "${path_current}data/train_data_aa.fa"

cp "${path_aa}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_aa.pkl"  # 不一定是terms_all,
cp "${path_aa}data/predictions.pkl" "${path_current}data/predictions_aa.pkl"
cp "${path_aa}step0_TrainTestSetting_local.py" "${path_current}step0_TrainTestSetting_local_aa.py"
cp "${path_aa}data/test_data.pkl" "${path_current}data/test_data_aa.pkl"
cp "${path_aa}data/test_data.fa" "${path_current}data/test_data_aa.fa"
cp "${path_aa}data/model_checkpoint.pth" "${path_current}data/model_checkpoint_aa.pth"

# ss8
#cd "$path_ss8" || exit 1
#python ../../s2_TrainTest/step4_Diamond.py  # 有一个aa的 diamond_aa.res 即可
cp "${path_ss8}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_ss8.pkl"  # 不一定是terms_all
cp "${path_ss8}data/predictions.pkl" "${path_current}data/predictions_ss8.pkl"
cp "${path_ss8}step0_TrainTestSetting_local.py" "${path_current}step0_TrainTestSetting_local_ss8.py"
cp "${path_ss8}data/test_data.pkl" "${path_current}data/test_data_ss8.pkl"
cp "${path_ss8}data/test_data.fa" "${path_current}data/test_data_ss8.fa"
cp "${path_ss8}data/model_checkpoint.pth" "${path_current}data/model_checkpoint_ss8.pth"


# 结合二者 (predictions_aa.pkl + predictions_ss8.pkl = predictions_aa_ss8.pkl)
#cd $path_current || exit 1
python AlphaBeta_Combine_aass8_PredictionsFile.py


