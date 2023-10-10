##!/bin/bash


# 当前文件夹要重新命名，i.e. step1_FindAlphaBeta_TrainMYCTU_TestMYCTU_aaK64F64_ss8K128F128


# 指定文件夹, aa & ss8
mkdir -p data
path_current="$(pwd)/"  # path_current后面没有/

path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"
path_aa="${path_base}output/test_TrainMYCTU_TestMYCTU_aa_DeepSS2GO_Kernel64_Filter64_Ontsall/"
path_ss8="${path_base}output/test_TrainMYCTU_TestMYCTU_ss8_DeepSS2GO_Kernel128_Filter128_Ontsall/"


# Diamond & 预处理文件
# aa
cd "$path_aa" || exit 1
python ../../s2_TrainTest/step4_Diamond.py  # 要在 path_aa中进行step4_Diamond.py
cp "${path_aa}data/diamond_aa.res" "${path_current}data/diamond_aa.res"
cp "${path_aa}data/train_data.pkl" "${path_current}data/train_data.pkl"


cp "${path_aa}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_aa.pkl"  # 不一定是terms_all,
cp "${path_aa}data/predictions.pkl" "${path_current}data/predictions_aa.pkl"

# ss8
#cd "$path_ss8" || exit 1
#python ../../s2_TrainTest/step4_Diamond.py  # 有一个aa的 diamond_aa.res 即可
cp "${path_ss8}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_ss8.pkl"  # 不一定是terms_all
cp "${path_ss8}data/predictions.pkl" "${path_current}data/predictions_ss8.pkl"


# 结合二者 (predictions_aa.pkl + predictions_ss8.pkl = predictions_aa_ss8.pkl)
cd $path_current || exit 1
python AlphaBeta_Combine_aass8_PredictionsFile.py



##########################################
####### Find AlphaBeta ###################
##########################################

## 这三个不一定全要做，根据具体情况来
python AlphaBeta_FindAlphasBetas.py -o bp
python AlphaBeta_FindAlphasBetas.py -o cc
python AlphaBeta_FindAlphasBetas.py -o mf


##########################################
########### evaluate #####################
##########################################




##########################################
########### predict #####################
##########################################




