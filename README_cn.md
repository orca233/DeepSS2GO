# DeepSS2GO_AB


## Cross species
跨物种训练了6类：ARATH, ECOLI, HUMAN, MOUSE, MYCTU, YEAST, & ALL00  
s1_DataPreprocessing_CrossSpecies/: 预处理数据  

s2_TrainTest/:   
- step1_SplitTrainTest_Terms: 拆分原始数据为train & test，找到对应terms 
- step2_Train: 训练，得到model + training.csv  
- step3_Test: 测试，得到predictions.pkl  
- step4_pkl2fa: 转格式
- step5_Diamond4CrossSpecies: 得到diamond_aa.res
- step6_FindAlpha: 找到最佳Alpha
- step7_EvaluateAlpha: 根据上面得到的Alpha，评估三大指标，Fmax, AUPR, Smin
- step8_PredictAlpha: 预测 test_data.pkl的结果，得到 results_bp/cc/mf.csv

>预测新的未知的数据有两种情况，
> 1. 只用aa/ss8的model。 (下面这两行命令用于情况1)  
> 2. 用aa+ss8的model。 (详见章节 PredictNew_AlphaBeta)  

step9_Diamond4New: 为全新的，未知的数据 更新Diamond  
step10_Predict_New: 根据aa/ss8 对未知数据进行预测  




## CAFA3
可以横向和其他文章做对比。


## TrainALL & Predict New

前两步**相同**：   
1. 准备初始fasta文件：pub_data/data_new/new_aa.fa  
2. 执行 s1_DataPreprocessing_New中的step1-7，得到在pub_data/data_new/中：  
new_clean_aa.pkl  
new_clean_ss8.pkl

### 只用 aa/ss8_model + Diamond 预测

3. 把这两个pkl文件cp到对应文件夹中。  
e.g. output/test_TrainALL00_TestALL00_aa_DeepSS2GO_Kernel8_Filter65536_Ontsall/data/

4. 执行：  
step9_Diamond4New: 为全新的，未知的数据 更新Diamond  
step10_Predict_New: 根据aa/ss8 对未知数据进行预测


### 用 aa_model + ss8_model + Diamond 预测
3. 数据预处理  
cd PredictNew_AlphaBeta/s1_DataPreprocessing_PredictNew/  

执行该文件夹下:  
- step1_fa2pkl: 转格式 
- step2_New_x_SPOT1DLM: 











