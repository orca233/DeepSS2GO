## 如果是个新文件，未知的，待测的，要测GO。
## 1）执行 s1_DataPreprocessing_new的 step1-8，得到 new_clean_aa/ss8.fa/pkl
## 2）拷贝到本文件夹的 data/
## 3）预测，step9 或 step10运行一个即可，因为一个文件夹下只有一个用aa或ss8训练的model.pth
#
##path_current="$(pwd)/"  # path_current后面没有/
#path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"
#
#
#cp "${path_base}pub_data/data_new/new_clean_aa.fa" data/
#cp "${path_base}pub_data/data_new/new_clean_aa.pkl" data/
#cp "${path_base}pub_data/data_new/new_clean_ss8.fa" data/
#cp "${path_base}pub_data/data_new/new_clean_ss8.pkl" data/
#
#
#
## diamond
###### 需要重新跑Diamond，因为这个是用 data_new/new_clean_aa.fa作为查询querying的 ！！！！！！！！！！！！！！
## 在当下文件夹 Diamond，因为 不是用 aa/ 文件夹中的test.fa，而是用这个New.fa
## Diamond
#echo ......... a long, long time ago ..........
#echo starting Diamond for new/unknown data
#
## 建库， train_data.fa--> train_data.dmnd 库
#echo -------------- build database train_data.dmnd --------------
#diamond makedb --in data/train_data.fa -d data/train_data
#
## 比对 查询 querying
#echo -------------- querying from database, creating diamond_aa.res --------------
##diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # routine
#diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/new_clean_aa.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # New





# 预测
echo predicting bp
python AlphaBeta_PredictAlphaBeta.py -t 0.02 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o bp  # 0.6

#echo predicting cc
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o cc  # 0.6
#
#echo predicting mf
#python AlphaBeta_PredictAlphaBeta.py -t 0.1 --in-file-aa 'data/new_clean_aa.fa' --in-file-ss8 'data/new_clean_ss8.fa' --test-data-file-aa 'data/new_clean_aa.pkl' --test-data-file-ss8 'data/new_clean_ss8.pkl' --alpha 'json' --beta 'json' -o mf  # 0.6







