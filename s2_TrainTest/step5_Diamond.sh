
#第二步，Diamond  (不论 aa / ss，都是用 aa 做diamond比对)，因为diamond的种子算法只认20个氨基酸，不认8个或3个的ss

#-d          = train_data.dmnd (数据库)
#-q (query)  = test_data.fa （要查询的文件）
#-o (output) = test_diamond.res （输出文件）
#qseqid      = Query Seq-id    查询列
#qseqid      = Subject Seq-id  数据库列
#bitscore    = Bit score     分数

#IN:
#train_data.dmnd    Diamond数据库
#test_data.fa            查询文件

#OUT:
#diamond_aa.res


echo ......... a long, long time ago ..........
echo starting step5_Diamond

# 建库， train_data.fa--> train_data.dmnd 库
echo -------------- build database train_data.dmnd --------------
diamond makedb --in data/train_data.fa -d data/train_data

# 比对 查询 querying
echo -------------- querying from database, creating diamond_aa.res --------------
diamond blastp -d data/train_data.dmnd --more-sensitive -t /tmp -q data/test_data.fa --outfmt 6 qseqid sseqid bitscore -o data/diamond_aa.res



echo .......... And they all lived happily ever after! ..........



