# 如果是个新文件，未知的，待测的，要测GO。


# 2）if 用model_aa预测，转pkl，准备好new_aa.fa & new_aa.pkl即可predict
# 3）if 用model_ss8预测，运行？？？？？？？
path_current="$(pwd)/"  # path_current后面没有/


# 1）首先给出 new_aa.fa，存在 data/中

# 2) fa2pkl: new_aa.pkl中有三列：['index', 'sequences', 'prop_annotations']
python Alpha_fa2pkl.py -inf data/new_aa.fa -ouf data/new_aa.pkl

# 2) s1_DataPreprocessing_New/中一步一步过
path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"

cd ${path_base} || exit
mkdir -p data_new

cd ${path_current} || exit
cp data/new_aa.fa "${path_base}pub_data/data_new/"
cp data/new_aa.pkl "${path_base}pub_data/data_new/"


# 接下来要在该文件夹下执行一系列的py
cd "${path_base}s1_DataPreprocessing_New/" || exit
#python step2...py




#echo predicting bp
#python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o bp --alpha 'json' # 0.6
#
#echo predicting cc
#python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o cc --alpha 'json' # 0.6
#
#echo predicting mf
#python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_aa.fa' --test-data-file 'data/new_aa.pkl' -o mf --alpha 'json' # 0.6







