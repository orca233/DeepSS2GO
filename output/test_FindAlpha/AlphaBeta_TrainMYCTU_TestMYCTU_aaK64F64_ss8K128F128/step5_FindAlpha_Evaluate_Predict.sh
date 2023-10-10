
mkdir -p results
mkdir -p results_alpha
mkdir -p results_alphabeta


#####################################################
# find_alpha  只是找到最佳alpha
# 根据实际情况，不一定要计算所有的 ont

echo Find Alpha bp
python Alpha_FindAlpha.py -o bp -ncc 60 --alpha-range 25,80,1  # -ar 三个数字之间不要有空格
echo Find Alpha cc
python Alpha_FindAlpha.py -o cc -ncc 60 --alpha-range 25,80,1
echo Find Alpha mf
python Alpha_FindAlpha.py -o mf -ncc 60 --alpha-range 25,80,1


#####################################################
# 前面 train & test 中，已经从model的到了 predictions.pkl，在test_data.pkl 后加上了两列，
# labels & preds [0, 0, 1, 0.....], [0.0045, 0.25, 0.001, ...]

# Evaluate, 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 从predictions.pkl里分析 precision - recall 曲线

# 如果alpha=NA，则采用json数据，否则alpha=数字，或外来click引入，
# alpha 为deepgo预测的
# alpha=0，全由diamond主导， alpha=1，全为DeepSS2GO深度学习主导

echo evaluate bp
python Alpha_Evaluate.py -o bp -a json
echo evaluate cc
python Alpha_Evaluate.py -o cc -a json
echo evaluate mf
python Alpha_Evaluate.py -o mf -a json

#####################################################
# Predict ， 当alpha=1 时，即全有deepSS2GO统计，alpha=0时，即全由diamond统计
# 需要 model_checkpoint.pth

echo predict bp
python Alpha_Predict.py -t 0.1   # t=threshold
















