
echo ......... a long, long time ago ..........
echo starting step5_Diamond

# build lib
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data.fa -d data/train_data

# querying
echo -------------- querying from database, creating diamond_aa.res --------------
diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res



echo .......... And they all lived happily ever after! ..........



