##!/bin/bash

# Diamond
echo ......... a long, long time ago ..........
echo starting Diamond for new/unknown data

# train_data.fa--> train_data.dmnd
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data_aa.fa -d data/train_data_aa

#  querying
echo -------------- querying from database, creating diamond_aa.res --------------
diamond blastp -d data/train_data_aa.dmnd --more-sensitive -t /tmp -q data/test_data_aa.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res




echo .......... And they all lived happily ever after! ..........





