echo ......... a long, long time ago ..........
echo starting Diamond for new/unknown data

# 建库， train_data.fa--> train_data.dmnd 库
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data.fa -d data/train_data

# 比对 查询 querying
echo -------------- querying from database, creating diamond_aa.res --------------
#diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # routine
diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # New
