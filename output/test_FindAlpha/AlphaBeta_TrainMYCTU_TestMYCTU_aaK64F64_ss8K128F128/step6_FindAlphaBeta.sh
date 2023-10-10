##!/bin/bash


# 新建一个文件夹，用来存放更新数据。下面的命令是在这个已经建好的文件夹中进行的。
# 1）e.g. mkdir AlphaBeta_TrainMYCTU_TestMYCTU_aaK64F64_ss8K128F128
# 2) cp Alpha*.py & step6_FindAlphaBeta.sh 到当前文件夹
# 3) 修改path_aa & path_ss8



# 指定文件夹, aa & ss8
mkdir -p data
path_current="$(pwd)/"  # path_current后面没有/

path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"
path_aa="${path_base}output/test_TrainMYCTU_TestMYCTU_aa_DeepSS2GO_Kernel64_Filter64_Ontsall/"
path_ss8="${path_base}output/test_TrainMYCTU_TestMYCTU_ss8_DeepSS2GO_Kernel128_Filter128_Ontsall/"



# Diamond & 预处理文件
# aa
cd "$path_aa" || exit 1
#python step4_Diamond.py
cp "${path_aa}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_aa.pkl"  # 不一定是terms_all,
cp "${path_aa}data/predictions.pkl" "${path_current}data/predictions_aa.pkl"

# ss8
cd "$path_ss8" || exit 1
#python step4_Diamond.py
cp "${path_ss8}data/terms_gominre_trxte.pkl" "${path_current}data/terms_gominre_trxte_ss8.pkl"  # 不一定是terms_all
cp "${path_ss8}data/predictions.pkl" "${path_current}data/predictions_ss8.pkl"


cd $path_current || exit 1
python AlphaBeta_Combine_aass8_PredictionsFile.py

