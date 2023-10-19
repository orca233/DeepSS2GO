


#####################################################
# Predict ， 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 需要 model_checkpoint.pth

# 这应该是先有的fa，然后转化出了

# 指定alpha是否=0/1 怎么改？？？？


# in-file & test-data-file 是可以更换成"待测的未知的全新的"
# 每一次只选bp/cc/mf中的一个来预测，因为其他两个默认为0




echo predicting bp
python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o bp  # 0.6

#echo predicting cc
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o cc  # 0.6
#
#echo predicting mf
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/test_data_aa.fa' --in-file-ss8 'data/test_data_ss8.fa' --test-data-file-aa 'data/test_data_aa.pkl' --test-data-file-ss8 'data/test_data_ss8.pkl' --alpha 'json' --beta 'json' -o mf  # 0.6








