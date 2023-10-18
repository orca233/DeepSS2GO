# DeepSS2GO_AB


## Train & Test




## CAFA3





## Predict New
不论预测 aa 还是 ss8，
1）准备初始fasta文件：pub_data/data_new/new_aa.fa

2）执行 s1_DataPreprocessing_New中的step1-7，得到在pub_data/data_new/中：
new_clean_aa.pkl
new_clean_ss8.pkl


把这两个cp到带model.pth的文件夹中。执行Alpha_PredictAlpha.py


