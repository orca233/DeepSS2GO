# DeepSS2GO
三大类：   
CrossSpecies/ 不同物种之间的交叉训练验证     
CAFA3/ 横向对比三大指标结果  
PredictNew/ 给出新的未知的fasta，GO预测  
其他文件夹：  
output/ 存放CS中大量的输出文件  
pub_data/ 存放三大类的公共数据  
redundancy/ 存放 data preprocessing中生成的大量 *.npt & *.csv 文件  


## Cross species
跨物种训练6类：ARATH, ECOLI, HUMAN, MOUSE, MYCTU, YEAST, & ALL00  

![CrossSpecies_Stage1_DataPreprocessing](figs/CrossSpecies_Stage1_DataPreprocessing.png)

s1_DataPreprocessing_CrossSpecies/: 预处理数据  
- step1_Uniprot2Swissprot.py: 挑出手工注释蛋白质
- step2_Swissprot_x_SPOT1DLM.py: 删除长度超过1024的蛋白质，和非标准氨基酸
- step3_SPOT1DLM_generage_esm.py: 生成 *esm.npy
- step4_SPOT1DLM_generate_prottrans.py: 生成 *pt.npy
- step5_SPOT1DLM_run_inference.py： 生成 *.csv
- step6_SPOT1DLM_csv_2_aass3ss8.py: 生成 SPOT1DLM_aass3ss8.pkl 
- step7_aa_2_ss3ss8.py: 生成 swissprot_clean_ALL00_ss3.pkl & swissprot_clean_ALL00_ss8.pkl
- step8_ClassifySpecies.py: 生成不同物种的 ss3/ss8.pkl


![CrossSpecies_Stage1_DataPreprocessing](figs/CrossSpecies_Stage2_TrainTest.png)

s2_TrainTest/: 训练，测试，评估   
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

- step9_Diamond4New: 为全新的，未知的数据 更新Diamond  
- step10_Predict_New: 根据aa/ss8 对未知数据进行预测  

> s2_TrainTest 使用方法：

CrossSpecies/s2_TrainTest/ 为global文件夹，  

在构建好pub_data后，比如需要 TrainECOLI_TestYEAST_aa，则
1. 在 CrossSpecies/下将整个 s2_TrainTest/ 复制为：
   s2_TrainTest_TrainECOLI_TestYEAST_aa/ (此为global文件夹)
2. 修改新文件夹中的 step0_TrainTestSetting_global.py 参数。  
   特别注意： 
   - 'device_ids': [0, 1]
   - 'dir0': 'test_TrainHUMAN_TestHUMAN_aa_test/'
   - 'aa_ss': 'aa'
   - 'train_data': 'HUMAN'
   - 'test_data': 'HUMAN'
   - 'kernels': [8, 16, 24, ...]
   - 'filters': [16, 32, 64, ...]   
   
3. 执行 run_KernelX_FilterY.py，会在 output/中新建文件夹：  
   test_TrainECOLI_TestYEAST_aa/ (此为local文件夹的母文件夹 = dir0)  
   并结合不同的kernel&filter生成一系列的local子文件夹：  
   e.g. DeepSS2GO_Kernel104_Filter8192_Ontsall/ (此为local文件夹的子文件夹 = dir1)  
   并将global文件夹中的文件cp到local子文件夹中，生成新的step0_TrainTestSetting_local.py  
   i.e.  
   cp -rf CrossSpecies/s2_TrainTest_TrainX_TestY_aa/ output/dir0/dir1/  
   e.g.  
   cp -rf CrossSpecies/s2_TrainTest_TrainECOLI_TestYEAST_aa/ output/test_TrainECOLI_TestYEAST_aa/DeepSS2GO_Kernel104_Filter8192_Ontsall/  
4. 在每一个local子文件夹中，根据step0_TrainTestSetting_local.py的具体参数，执行: step1-8。
   一般可能只运行 step1-3，后面选择合适的单独运行 step4-8用来 find_alpha。


## CAFA3
可以横向和其他文章做对比。







## TrainALL & Predict New

### case 0: 只用 aa/ss8_model + Diamond 预测
前两步**相同**：   
1. 准备初始fasta文件，aa格式：pub_data/data_new/new_aa.fa  
2. 执行 s1_DataPreprocessing_New中的step1-7，得到在pub_data/data_new/中：  
new_clean_aa.pkl  
new_clean_ss8.pkl







3. 把这两个pkl文件cp到对应文件夹中。  
e.g. output/test_TrainALL00_TestALL00_aa_DeepSS2GO_Kernel8_Filter65536_Ontsall/data/

4. 执行：  
step9_Diamond4New: 为全新的，未知的数据 更新Diamond  
step10_Predict_New: 根据aa/ss8 对未知数据进行预测


### case 1: 用 aa_model + ss8_model + Diamond 预测
3. 数据预处理  
cd PredictNew_AlphaBeta/s1_DataPreprocessing_PredictNew/  

执行该文件夹下:  
- step1_fa2pkl: 转格式 
- step2_New_x_SPOT1DLM: 















































