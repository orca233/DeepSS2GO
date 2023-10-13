# 如果是个新文件，未知的，待测的，要测GO。
# 1）首先给出 new_aa.fa，存在 data/中
# 2）if 用model_aa预测，转pkl，准备好new_aa.fa & new_aa.pkl即可predict
# 3）if 用model_ss8预测，运行？？？？？？？




python Alpha_fa2pkl.py -inf data/new_aa.fa -ouf data/new_aa.pkl


echo predicting bp
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o bp --alpha 'json' # 0.6

echo predicting cc
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o cc --alpha 'json' # 0.6

echo predicting mf
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o mf --alpha 'json' # 0.6







