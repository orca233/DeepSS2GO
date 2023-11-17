# 如果是个新文件，未知的，待测的，要测GO。
# 1）执行 s1_DataPreprocessing_new的 step1-8，得到 new_clean_aa/ss8.fa/pkl
# 2）拷贝到本文件夹的 data/
# 3）预测，step9 或 step10运行一个即可，因为一个文件夹下只有一个用aa或ss8训练的model.pth

path_current="$(pwd)/"  # path_current后面没有/
path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO_AB/"
#path_base=$(cd ../../ && pwd)
#path_base="$(path_base)/"

cp "${path_base}pub_data/data_new/new_clean_aa.fa" data/
cp "${path_base}pub_data/data_new/new_clean_aa.pkl" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.fa" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.pkl" data/



# diamond
##### 需要重新跑Diamond，因为这个是用 data_new/new_clean_aa.fa作为查询querying的 ！！！！！！！！！！！！！！
# 在当下文件夹 Diamond，因为 不是用 aa/ 文件夹中的test.fa，而是用这个New.fa
# Diamond
echo ......... a long, long time ago ..........
echo starting Diamond for new/unknown data

# 建库， train_data.fa--> train_data.dmnd 库
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data.fa -d data/train_data

# 比对 查询 querying
echo -------------- querying from database, creating diamond_aa.res --------------
#diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # routine
diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/new_clean_aa.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # New


#这里要查看输出结果，有多少匹配上了!!!!!!!!!!!!!!!!!!!!!
#e.g.
#Reported 45470 pairwise alignments, 45470 HSPs.
#3252 queries aligned.



