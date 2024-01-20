##!/bin/bash


# 在当下文件夹 Diamond，因为 不是用 aa/ 文件夹中的test.fa，而是用这个New.fa
# Diamond
echo ......... a long, long time ago ..........
echo starting Diamond for new/unknown data

# 建库， train_data.fa--> train_data.dmnd 库
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data_aa.fa -d data/train_data_aa

# 比对 查询 querying
echo -------------- querying from database, creating diamond_aa.res --------------
diamond blastp -d data/train_data_aa.dmnd --more-sensitive -t /tmp -q data/test_data_aa.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res  # original, routine




echo .......... And they all lived happily ever after! ..........





