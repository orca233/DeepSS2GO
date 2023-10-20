##!/bin/bash


# 当前文件夹要重新命名，i.e. step1_FindAlphaBeta_TrainMYCTU_TestMYCTU_aaK64F64_ss8K128F128



# 结合二者 (predictions_aa.pkl + predictions_ss8.pkl = predictions_aa_ss8.pkl)
#cd $path_current || exit 1
python AlphaBeta_Combine_aass8_PredictionsFile.py



##########################################
####### Find AlphaBeta ###################
##########################################

## 这三个不一定全要做，根据具体情况来。分两轮，这个 '0, 101, 10'在第二轮要修改
#python AlphaBeta_FindAlphaBeta.py -ncc 80 -o bp --alpha-range '0, 101, 10' --beta-range '0, 101, 10'
#python AlphaBeta_FindAlphaBeta.py -ncc 80 -o cc --alpha-range '0, 101, 10' --beta-range '0, 101, 10'
#python AlphaBeta_FindAlphaBeta.py -ncc 80 -o mf --alpha-range '0, 101, 10' --beta-range '0, 101, 10'



########### 第二轮 #####################
# 因为 bp 在a=0.2, b=0.3取最大，所以新的range 去 0.1-0.3, 0.2-0.4，间隔均为2：

python AlphaBeta_FindAlphaBeta.py -ncc 80 -o bp --alpha-range '10, 31, 2' --beta-range '20, 41, 2'
#python AlphaBeta_FindAlphaBeta.py -ncc 80 -o cc --alpha-range '0, 101, 2' --beta-range '0, 101, 2'
#python AlphaBeta_FindAlphaBeta.py -ncc 80 -o mf --alpha-range '0, 101, 2' --beta-range '0, 101, 2'


# 得到 (0.16, 0.28, -0.535361639858736)




