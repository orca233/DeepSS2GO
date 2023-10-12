

#####################################################
# Predict ， 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 需要 model_checkpoint.pth

echo predict bp
python Alpha_Predict.py -t 0.1   # t=threshold  应该可以指定alpha是否=0/1 怎么改？？？？？






