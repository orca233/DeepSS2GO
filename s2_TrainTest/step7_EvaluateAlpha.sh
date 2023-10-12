

#####################################################
# 前面 train & test 中，已经从model的到了 predictions.pkl，在test_data.pkl 后加上了两列，
# labels & preds [0, 0, 1, 0.....], [0.0045, 0.25, 0.001, ...]

# Evaluate, 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 从predictions.pkl里分析 precision - recall 曲线

# 如果alpha=json，则采用json数据，否则alpha=数字，或外来click引入，
# alpha 为deepgo预测的
# alpha=0，全由diamond主导， alpha=1，全为DeepSS2GO深度学习主导，alpha可设为0.3等float

echo evaluate bp
python Alpha_Evaluate.py -o bp --alpha 'json' # 0.6
#echo evaluate cc
#python Alpha_Evaluate.py -o cc --alpha json
#echo evaluate mf
#python Alpha_Evaluate.py -o mf --alpha json





