

#####################################################
# Predict ， 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 需要 model_checkpoint.pth

# 这应该是先有的fa，然后转化出了

# 指定alpha是否=0/1 怎么改？？？？


# in-file & test-data-file 是可以更换成"待测的未知的全新的"

echo predicting
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/test_data.fa' --test-data-file 'data/test_data.pkl'






