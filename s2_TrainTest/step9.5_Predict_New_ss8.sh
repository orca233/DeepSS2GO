# 如果是个新文件，未知的，待测的，要测GO。
# 1）执行 s1_DataPreprocessing_new的 step1-8，得到 new_clean_aa/ss8.fa/pkl
# 2）拷贝到本文件夹的 data/
# 3）预测，step9 或 step10运行一个即可，因为一个文件夹下只有一个用aa或ss8训练的model.pth

#path_current="$(pwd)/"  # path_current后面没有/
path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"


cp "${path_base}pub_data/data_new/new_clean_aa.fa" data/
cp "${path_base}pub_data/data_new/new_clean_aa.pkl" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.fa" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.pkl" data/


# 预测
echo predicting bp
python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o bp --alpha 'json' # 0.6
#
#echo predicting cc
#python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o cc --alpha 'json' # 0.6
#
#echo predicting mf
#python Alpha_PredictAlpha.py -t 0.1 --in-file 'data/new_clean_ss8.fa' --test-data-file 'data/new_clean_ss8.pkl' -o mf --alpha 'json' # 0.6







