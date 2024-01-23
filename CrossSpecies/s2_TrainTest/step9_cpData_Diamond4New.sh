

path_current="$(pwd)/"
path_base="/home/fsong/work/py_proj/prot_algo/DeepSS2GO/"


cp "${path_base}pub_data/data_new/new_clean_aa.fa" data/
cp "${path_base}pub_data/data_new/new_clean_aa.pkl" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.fa" data/
cp "${path_base}pub_data/data_new/new_clean_ss8.pkl" data/




# Diamond
echo ......... a long, long time ago ..........
echo starting Diamond for new/unknown data

# train_data.fa--> train_data.dmnd
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data.fa -d data/train_data

# querying
echo -------------- querying from database, creating diamond_aa.res --------------
diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/new_clean_aa.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # New





