

#####################################################
# Predict ， 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 需要 model_checkpoint.pth

# 这应该是先有的fa，然后转化出了

# 指定alpha是否=0/1 怎么改？？？？


# in-file & test-data-file 是可以更换成"待测的未知的全新的"
# 每一次只选bp/cc/mf中的一个来预测，因为其他两个默认为0

# (1 - alpha - beta) * diamond + alpha * preds_aa + beta * preds_ss8
# alpha=0: 全由 diamond 统计
# alpha=1: 全由 deepSS2GO_aa 统计
# beta=1: 全由 deepSS2GO_ss8 统计

echo predicting bp
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o bp --alpha 'json' # 0.6

echo predicting cc
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o cc --alpha 'json' # 0.6

echo predicting mf
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl' -o mf --alpha 'json' # 0.6






