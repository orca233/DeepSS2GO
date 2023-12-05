# PredictNew

[TOC]



## 数据预处理，将aa.fa转成ss8.fa并生成*.pkl

进入主目录：/PredictNew/s1_DataPreprocessing_PredictNew/

### 1. 下载预训练模型(step3/4用到)：
- esm1b_t33_650M_UR50S  # 自动下载超慢
- Prot_T5_XL_UniRef50

下载链接：
```bash
https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
# 保存在：/home/USERNAME/.cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S-contact-regression.pt & esm1b_t33_650M_UR50S.pt

https://huggingface.co/Rostlab/prot_t5_xl_uniref50/tree/main
# 保存在指定目录：/home/USERNAME/.../Prot_T5_XL_UniRef50/

# 修改step0_DataPreprocessingSetting.py中
path_Prot_T5_XL_UniRef50 = /home/USERNAME/.../Prot_T5_XL_UniRef50/

```


### 2. 准备初始fasta文件，
aa格式，存放在：/pub_data/data_new/new_aa.fa

比如：
    
```
>slam1
MVIFYFCGKTFMPARNRWMLLLPLLASAAYAEETPREPDLRSRPEFRLHEAEVKPIDREKVPGQVREKGKVLQIDGETLLKNPELLSRAMYSAVVSNNIAGIRVILPIYLQQAQQDKMLALYAQGILAQADGRVKEAISHYRELIAAQPDAPAVRMRLAAALFENRQNEAAADQFDRLKAENLPPQLMEQVELYRKALRERDAWKVNGGFSVTREHNINQAPKRQQYGKWTFPKQVDGTAVNYRLGAEKKWSLKNGWYTTAGGDVSGRVYPGNKKFNDMTAGVSGGIGFADRRKDAGLAVFHERRTYGNDAYSYTNGARLYFNRWQTPKWQTLSSAEWGRLKNTRRARSDNTHLQISNSLVFYRNARQYWMGGLDFYRERNPADRGDNFNRYGLRFAWGQEWGGSGLSSLLRLGAAKRHYEKPGFFSGFKGERRRDKELNTSLSLWHRALHFKGITPRLTLSHRETRSNDVFNEYEKNRAFVEFNKTF
>slam2
MLYFRYGFLVVWCAAGVSAAYGADAPAILDDKALLQVQRSVSDKWAESDWKVENDAPRVVDGDFLLAHPKMLEHSLRDALNGNQADLIASLADLYAKLPDYDAVLYGRARALLAKLAGRPAEAVARYRELHGENAADERILLDLAAAEFDDFRLKSAERHFAEAAKLDLPAPVLENVGRFRKKTEGLTGWRFSGGISPAVNRNANNAAPQYCRQNGGRQICSVSRAERAAGLNYEIEAEKLTPLADNHYLLFRSNIGGTSYYFSKKSAYDDGFGRAYLGWQYKNARQTAGILPFYQVQLSGSDGFDAKTKRVNNRRLPPYMLAHGVGVQLSHTYRPNPGWQFSVALEHYRQRYREQDRAEYNNGRQDGFYVSSAKRLGESATVFGGWQFVRFVPKRETVGGAVNNAAYRRNGVYAGWAQEWRQLGGLNSRVSASYARRNYKGIAAFSTEAQRNREWNVSLALSHDKLSYKGIVPALNYRFGRTESNVPYAKRRNSEVFVSADWRF

```


### 3. 执行命令

执行 s1_DataPreprocessing_New中的step1-8，最终在pub_data/data_new/生成：
- new_clean_aa.pkl 
- new_clean_aa.fa
- new_clean_ss8.pkl 
- new_clean_ss8.fa


## 单独预测 aa or ss8，仅用到 alpha

> 前情概要，在结束Analysis_ALL00_CAFA3_FindAlpha_FindAlphaBeta.md后，我们得到了两个文件夹：
- /Analysis/ALL00_CAFA3_FindAlpha_aaORss8/
- /Analysis/ALL00_CAFA3_FindAlphaBeta_aaANDss8/


### 1. 同步文件夹
把 /Analysis/ALL00_CAFA3_FindAlpha_aaORss8/ 的三个ALL00文件夹
rsync到：/PredictNew/s2_PredictNew_Alpha/
> ALL00_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall/
> ALL00_ss8_DeepSS2GO_Kernel32_Filter32768_Ontsall/
> ALL00_ss8_DeepSS2GO_Kernel48_Filter16384_Ontsall/


### 2. 预测

以 ALL00_aa_DeepSS2GO_Kernel16_Filter32768_Ontsall/ 为例：
> 同理执行 
> ALL00_ss8_DeepSS2GO_Kernel32_Filter32768_Ontsall/
> ALL00_ss8_DeepSS2GO_Kernel48_Filter16384_Ontsall/

```bash
# 修改 step9_cpData_Diamond4New.sh 中为个人账户路径
path_base="/home/USERNAME/work/py_proj/prot_algo/DeepSS2GO/"

bash step9_cpData_Diamond4New.sh  # 把这4个pkl/fa文件cp到对应文件夹中，并diamond

python step10_Predict_New.py  # threshold 改成0.02 !!!!!!

```

输出的结果在相应的 /data/里：
results_bp/cc/mf.csv


<font color=red size="5">**特别注意：**</font>   

ALL00_aa/ 里的 results_bp/cc/mf皆可用
ALL00_ss8_K32_F32768/ 只看 results_bp/mf
ALL00_ss8_K48_F16384/ 只看 results_cc


### 3. 看结果

以slam1/2为例，看看结果

对于bp & mf，左侧是aa预测，右侧是ss8的预测。
可以看到 ss8 预测里面有transport

<font color=green size="5">**有意思的是：**</font>

slam1 & 2 在Diamond中是：0 queries aligned.


![Img](/Users/sujiaqi/Desktop/SF-readme/WechatIMG1340.jpg)

![Img](/Users/sujiaqi/Desktop/SF-readme/WechatIMG1341.jpg)





## 结合 aa & ss8，用到 alpha + beta


> 前情概要，在结束Analysis_ALL00_CAFA3_FindAlpha_FindAlphaBeta.md后，我们得到了两个文件夹：
- /Analysis/ALL00_CAFA3_FindAlpha_aaORss8/
- /Analysis/ALL00_CAFA3_FindAlphaBeta_aaANDss8/


### 1. 同步文件夹

把 /Analysis/ALL00_CAFA3_FindAlphaBeta_aaANDss8/ 的三个ALL00文件夹
rsync到：/PredictNew/s3_PredictNew_AlphaBeta/
> s3_AlphaBeta_TrainALL00_TestALL00_bp_aaK16F32768_ss8K32F32768/
> s3_AlphaBeta_TrainALL00_TestALL00_cc_aaK16F32768_ss8K48F16384/
> s3_AlphaBeta_TrainALL00_TestALL00_mf_aaK16F32768_ss8K32F32768/


### 2. 预测

以s3_AlphaBeta_TrainALL00_TestALL00_bp_aaK16F32768_ss8K32F32768/为例：

> 同理执行
> s3_AlphaBeta_TrainALL00_TestALL00_cc_aaK16F32768_ss8K48F16384/
> s3_AlphaBeta_TrainALL00_TestALL00_mf_aaK16F32768_ss8K32F32768/


```bash
# 修改 [step9_cpData_Diamond4New.sh](step6_cpData_Diamond4New.sh) 中为个人账户路径
path_base="/home/USERNAME/work/py_proj/prot_algo/DeepSS2GO/"

bash step6_cpData_Diamond4New.sh  # 把这4个pkl/fa文件cp到对应文件夹中，并diamond

bash step7_PredictAlphaBeta_New.sh  # threshold 改成0.02 !!!!!!

```


输出的结果在相应的 /data/里：
results_bp/cc/mf.csv


### 3. 看结果



#### bp

![Img](/Users/sujiaqi/Desktop/SF-readme/WechatIMG1342.jpg)



#### cc

![Img](/Users/sujiaqi/Desktop/SF-readme/WechatIMG1343.jpg)




#### mf

![Img](/Users/sujiaqi/Desktop/SF-readme/WechatIMG1344.jpg)







