## README 至少要保证所有的文件夹存在，否则输出的excel，该行会左移对其，就乱了
# 要保证kernel/filter 小的部分存在，否则csv乱码，会自动向左靠紧

# 有三部分，
# 第1部分，读取数据，生成 csv。
# 第2部分，plot
# 第3部分，count_epochs，如果epochs很高，说明lr设置的太小。如果epoch=7,8，说明lr太大了


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import os
import csv


#####################################################
########## 第1部分，读取数据，生成 csv #################
#####################################################


#################### pre_setting ################

# path_base = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'   # Keras 结果
path_base = '/home/fsong/work/py_proj/prot_algo/'  # Pytorch 结果

# 设置根目录
# dir0 = 'test_TrainHUMAN_TestHUMAN_ss8_MAXLEN1024/step3_done/'
# root_path = '%soutput/%s' % (path_base, dir0)

# 初始化kernel和filter的取值范围
# 2 ^ 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18

#################################################################
############################## pytorch ############################
#################################################################



############ CAFA3 ##################

# round_0
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3/test_CAFA3_aa/step3_done/'
# kernel_range = [4, 8, 12, 16, 24, 32]  # 其他的重新计算中。。。
# filter_range = [4096, 8192, 16384, 32768, 65536]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3/test_CAFA3_ss8/step3_done/'
# kernel_range = [24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]
# filter_range = [1024, 2048, 4096, 8192, 16384, 32768, 65536]


# round_1
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3_aa_reverse/step3_done/'
# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 104, 112, 120, 128, 256]  # 缺：64, 72, 80, 88, 96,
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3_ss8_reverse/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


# round_2  ---- 和round 1没啥区别，果然不应该运行 step1_split
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3_aa_reverse/round2_step3_done_lack_K64_72_80_88_96_just_recalc_step3/'
# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 104, 112, 120, 128, 256]  # 缺：64, 72, 80, 88, 96,
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


## CAFA3_round2 ------ 没有运行step1
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3_round_2_aa/'
# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3_round_2_ss8/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


## CAFA3_round3
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3/test_CAFA3_round3_aa_APUS/step3_done/'
# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_CAFA3/test_CAFA3_round3_ss8/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


# CAFA3_round5
# root_path = path_base + 'DeepSS2GO/output/test_CAFA3/test_CAFA3_round5_aa/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
#
# root_path = path_base + 'DeepSS2GO/output/test_CAFA3/test_CAFA3_round5_ss8/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


############## TrainALL00 ##################

### F=131072，由KernelX_X_Filter65536_65536通过GPU计算， 有K=4_4, 8_8, 12_12, 16_16, 24_24, 32_32, 40_40, 48_48，再大K算不了
### F=262144 由KernelX_X_Filter131072_131072通过CPU计算，有K=4_4, 8_8, 12_12, 16_16, 24_24

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainALL00/test_TrainALL00_TestALL00_aa/step3_done/'
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainALL00/test_TrainALL00_TestALL00_ss8/step3_done/'

# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144]

# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainALL00/test_TrainALL00_TestALL00_ss3/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # 没有 256
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]



# 新的ALL00 20231110
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00/test_TrainALL00_TestALL00_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00/test_TrainALL00_TestALL00_ss8/'


## FC = 1
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_aa_FC1/'
root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_ss8_FC1/'


# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
kernel_range = [8, 16, 24, 32, 40, 48]  # lack:
filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]





############### TrainARATH ################

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestARATH_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestARATH_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestECOLI_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestECOLI_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestHUMAN_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestHUMAN_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestMOUSE_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestMOUSE_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestMYCTU_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestMYCTU_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestYEAST_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainARATH_TestYEAST_ss8/'








############### TrainECOLI ################

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestARATH_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestARATH_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestECOLI_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestECOLI_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestHUMAN_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestHUMAN_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestMOUSE_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestMOUSE_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestMYCTU_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestMYCTU_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestYEAST_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainECOLI_TestYEAST_ss8/'



############### TrainHUMAN ################

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestARATH_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestARATH_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestECOLI_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestECOLI_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestHUMAN_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestHUMAN_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestMOUSE_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestMOUSE_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestMYCTU_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestMYCTU_ss8/'

# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestYEAST_aa/'
# root_path = path_base + 'DeepSS2GO/output/test_CrossSpecies/' + 'test_TrainHUMAN_TestYEAST_ss8/'


## FC = 1
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_aa_FC1/'
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_ss8_FC1/'

## FC = 2
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_aa_FC2/'
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_ss8_FC2/'



##########################



# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  # lack:
# # kernel_range = [8, 16, 24, 32]  # lack:
#
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]



# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestARATH_aa/'
#
# root_path = path_base + 'DeepSS2GO/output/test_TrainMOUSE_TestARATH_aa/'
#



# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestHUMAN_aa/step3_done/'
# kernel_range = [4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  # 有两个特殊点：K=4, 12
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144]
### F=262144 只有K=4,8,12,16,24，通过[131072, 131072]算出来的, cpu

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestHUMAN_ss8/step3_done/'  # 不需要 4 GPU算 F=131072，已经过拟合了
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072]

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestHUMAN_ss3/step3_done/'  # 不需要 4 GPU算 F=131072，已经过拟合了
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]



################################ 其他尝试 ###################################

# # HUMAN_HUMAN_aa_FCdepth1
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN_TestHUMAN_aa_MAXLEN1024_FCdepth1/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]  # F=262144 只有K=4,8,12,16,24

# # # # HUMAN_HUMAN_ss8_FCdepth1
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN_TestHUMAN_ss8_MAXLEN1024_FCdepth1/step3_done/'  # 还要 4 GPU算 F=131072
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

# # HUMAN_HUMAN_aa_MAXLEN2000
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN_TestHUMAN_aa_MAXLEN2000/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88]  # 有两个特殊点：K=4, 12
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]  # F=262144 只有K=4,8,12,16,24, F=65536双核跑不了

# # HUMAN_HUMAN_aa_MAXLEN2000
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN_TestHUMAN_aa_MAXLEN2000/step3_done/'
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80]  # 有两个特殊点：K=4, 12
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]  # F=262144 只有K=4,8,12,16,24, F=65536双核跑不了

################################ 其他尝试 end ###################################



################ TrainHUMAN_TestY ##################

# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestARATH_aa/step3_done/'
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestARATH_ss8/step3_done/'
# root_path = path_base + 'DeepSS2GO_Pytorch/output/test_TrainHUMAN/test_TrainHUMAN_TestARATH_ss3/step3_done/'


# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]












#################################################################
############################## Keras ############################
#################################################################



# # HUMAN_HUMAN_aa
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_KernelsXFiltersY_aa/'
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]


# # # HUMAN_HUMAN_ss8
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_KernelsXFiltersY_ss8/'  #
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]


# # # # HUMAN_HUMAN_aa_FCdepth1
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_KernelsXFiltersY_aa_FCdepth1/step3_done/'
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]


# # # # HUMAN_HUMAN_ss8_FCdepth1
# root_path = path_base + 'DeepSS2GO/output/test_TrainHUMAN_TestHUMAN_KernelsXFiltersY_ss8_FCdepth1/step3_done/'  #
# kernel_range = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]




# # ALL00_ALL00_aa
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_KernelsXFiltersY_aa/'
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]


# # ALL00_ALL00_ss8
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_KernelsXFiltersY_ss8/'  #
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]  #
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 49152]


# # # # ALL00_ALL00_aa_FCdepth1
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_KernelsXFiltersY_aa_FCdepth1/'
# kernel_range = [4, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512]
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]


# # # # # ALL00_ALL00_ss8_FCdepth1
# root_path = path_base + 'DeepSS2GO/output/test_TrainALL00_TestALL00_KernelsXFiltersY_ss8_FCdepth1/step3_K4_8_16/'  #
# kernel_range = [4, 8, 16]  # , 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256, 512
# filter_range = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384]  #, 32768, 49152


onts = 'all'  # onts 有几种形式： 'bp', 'cc', 'mf', 'all' 。 这个和 step0_global.py的统一


#################### start working ################


# 初始化表头
onts_list = []
if onts == 'all':
    onts_list = ['bp', 'cc', 'mf']
else:
    onts_list = [onts]

header = ['filter']
for kernel in kernel_range:
    for category in onts_list:  # ['bp', 'cc', 'mf'] 或 其中之一，为list形式
        header.extend([f'kernel{kernel}_{category}_Fmax',
                       f'kernel{kernel}_{category}_AUPR',
                       f'kernel{kernel}_{category}_Smin'])

# 初始化表格数据
table_data = []
for filter_size in filter_range:
    row_data = [filter_size]
    for kernel_size in kernel_range:
        folder_name = f"DeepSS2GO_Kernel{kernel_size}_Filter{filter_size}_Onts{onts}"
        # # 构建结果文件夹路径
        # if 'Pytorch' in root_path:
        #     folder_name = f"DeepSS2GO_Kernel{kernel_size}_Filter{filter_size}_Onts{onts}"  # 这里有3个txt  for Pytorch
        # else:
        #     folder_name = f"DeepSS2GO_Kernel{kernel_size}_Filter{filter_size}"  # 这里有3个txt   for Keras

        folder_path = os.path.join(root_path, folder_name, 'results_alpha')
        # 检查文件夹是否存在
        if not os.path.exists(folder_path):
            print(f"Folder {folder_path} does not exist")
            continue
        # 初始化每个kernel的6个值，即3个类别下的Fmax和AUPR
        fmax_bp, aupr_bp, smin_bp = '', '', ''
        fmax_cc, aupr_cc, smin_cc = '', '', ''
        fmax_mf, aupr_mf, smin_mf = '', '', ''




        # 遍历txt文件夹
        for filename in os.listdir(folder_path):
            # if filename.endswith('.txt'):
            if filename.endswith('1.00.txt'):
                # print(filename)

                # 提取类别
                category = filename.split('_')[3].split('.')[0]  # 不能用倒数第一个了，用第0123个 Fmax_AUPR_Smin_cc_1.00.txt  original: category = filename.split('_')[-1].split('.')[0]
                # 提取Fmax和AUPR信息
                with open(os.path.join(folder_path, filename), 'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if line.startswith('Fmax'):
                            fmax = line.split(':')[-1].strip()
                            if category == 'bp':
                                fmax_bp = fmax
                            elif category == 'cc':
                                fmax_cc = fmax
                            elif category == 'mf':
                                fmax_mf = fmax
                        elif line.startswith('AUPR'):
                            aupr = line.split(':')[-1].strip()
                            if category == 'bp':
                                aupr_bp = aupr
                            elif category == 'cc':
                                aupr_cc = aupr
                            elif category == 'mf':
                                aupr_mf = aupr

                        elif line.startswith('Smin'):
                            smin = line.split(':')[-1].strip()
                            if category == 'bp':
                                smin_bp = smin
                            elif category == 'cc':
                                smin_cc = smin
                            elif category == 'mf':
                                smin_mf = smin

        # 将结果添加到一行数据中


        if onts == 'all':
            # if filter_size == 262144:
                # print('2222222222222222')
                # print(folder_name)
                # print([fmax_bp, aupr_bp, smin_bp, fmax_cc, aupr_cc, smin_cc, fmax_mf, aupr_mf, smin_mf])

            row_data.extend([fmax_bp, aupr_bp, smin_bp, fmax_cc, aupr_cc, smin_cc, fmax_mf, aupr_mf, smin_mf])
            # print('row_data = ', row_data)

        elif onts == 'bp':
            row_data.extend([fmax_bp, aupr_bp, smin_bp])
        elif onts == 'cc':
            row_data.extend([fmax_cc, aupr_cc, smin_cc])
        elif onts == 'mf':
            row_data.extend([fmax_mf, aupr_mf, smin_mf])

    # 将一行数据添加到表格数据中
    table_data.append(row_data)
    # print('----- table_data = ', table_data, '\n')


# 将表头和表格数据写入csv文件
with open('CrossValidation_ColumnFilters.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)
    writer.writerows(table_data)



print('step1  done')



#####################################################
########## 第2部分，plot #############################
#####################################################

markers = ['|', '_', '.', ',', 'v', '^', '<', '>', 'x', 'o', 's', 'd', 'D', 'h', 'H', 'p', '1', '2', '3', '4', '+', '*']


# 读取csv文件
# data = pd.read_csv('CrossValidation_ColumnFilters2.csv')
data = pd.read_csv('CrossValidation_ColumnFilters.csv')


# 取第一列为横坐标
x = data.iloc[:, 0]
x = np.log2(x)

# 需要生成的九个数值对应的列索引
column_indices = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# column_names = ['bp_Fmax', 'bp_AUPR', 'bp_Smin', 'cc_Fmax', 'cc_AUPR', 'cc_Smin', 'mf_Fmax', 'mf_AUPR', 'mf_Smin']
column_names = ['bp_Fmax', 'cc_Fmax', 'mf_Fmax', 'bp_AUPR', 'cc_AUPR', 'mf_AUPR', 'bp_Smin', 'cc_Smin', 'mf_Smin']

# 生成一个大图，包含9个小图
fig, axs = plt.subplots(3, 3, figsize=(12, 12))

# Add a general title above the subplots  总标题
parts = root_path.split('/')
folder_name = parts[-2]  # 倒数第二个部分就是目标文件夹名称
print('folder_name=', folder_name)
fig.suptitle(folder_name, fontsize=16)


# 循环遍历每个小图的索引和对应的列名
for i, index in enumerate(column_indices):  # i 代表小图的索引
    row = i // 3  # 计算小图所在的行索引
    col = i % 3  # 计算小图所在的列索引

    # 在指定的子图上绘制数据
    ax = axs[row, col]

    # 初始化：每个小图内部的最大值和最小值
    min_value = 100
    max_value = 0
    max_x = 0  # 保存最大值对应的x坐标
    min_x = 0  # 保存最小值对应的x坐标


    # 取包含"bp_Fmax/cc_AUPR/mf_Smin"等的列，提取数据并作小图
    marker_index = 0
    for j in range(1, data.shape[1]):  # 逐一过csv中的每一列
        if column_names[i] in data.columns[j]:  # column_names[i]='bp_Fmax/cc_AUPR...'。 data.columns[j]是csv文件中第j列名字, eg:kernel8_cc_Smin
            # 获取小图对应的数据和标题
            y = data.iloc[:, j]
            title = data.columns[j]
            # print(marker_index)
            marker = markers[marker_index]  # 按顺序选择一个标记  markers[j % len(markers)]
            marker_index += 1
            ax.plot(x, y, label=title, linewidth=2, marker=marker, markersize=5)  # linestyle='-'


            # if i == 6 or i == 7 or i == 8:  # bp/cc/mf_smin 这个图，画legend
            # if 3 <= i and i <= 8:  # bp/cc/mf_smin 这个图，画legend
            ax.legend(fontsize=7, loc='lower left', ncol=2)  # 两列
            ax.set_xlabel('log2(Filter)', fontsize=10)
            ax.set_ylabel('Value', fontsize=10)
            ax.set_title(column_names[i], fontsize=12)  # 小图的标题
            # 设置横坐标轴的xticks间隔为1
            ax.set_xticks(np.arange(min(x) - 1, max(x) + 1, 1))




            # 更新最大值和最小值
            if np.max(y) > max_value:
                max_value = np.max(y)
                max_x = x[np.argmax(y)]
            if np.min(y) < min_value:
                min_value = np.min(y)
                min_x = x[np.argmin(y)]

    # 在小图内部添加文本显示最大值和最小值
    # 文本位置：
    text_x_max = max_x
    text_y_max = max_value
    text_x_min = min_x
    text_y_min = min_value

    # 获取小图的数据范围
    y_min, y_max = ax.get_ylim()
    x_min, x_max = ax.get_xlim()

    # 根据范围调整文本的位置
    text_y_max_offset = 0.02 * (y_max - y_min)
    text_y_min_offset = 0.02 * (y_max - y_min)
    text_x_offset = 0.1 * (x_max - x_min)


    # 在小图内部添加文本显示最大值和最小值
    ax.text(0.05, 0.95, f"Max: {max_value:.3f} (x={max_x:.1f})", fontsize=10, ha='left', va='top', transform=ax.transAxes)
    ax.text(0.05, 0.90, f"Min: {min_value:.3f} (x={min_x:.1f})", fontsize=10, ha='left', va='top', transform=ax.transAxes)

    # # 设置纵坐标轴的刻度间隔为Fmax&AUPR为0.01，对于Smin间隔为0.2
    # if column_names[i] == 'bp_Fmax':
    #     ax.set_yticks(np.arange(0.35, 0.43, 0.01))
    # if column_names[i] == 'cc_Fmax':
    #     ax.set_yticks(np.arange(0.58, 0.7, 0.01))
    # if column_names[i] == 'mf_Fmax':
    #     ax.set_yticks(np.arange(0.34, 0.48, 0.01))
    # if column_names[i] == 'bp_AUPR':
    #     ax.set_yticks(np.arange(0.28, 0.38, 0.01))
    # if column_names[i] == 'cc_AUPR':
    #     ax.set_yticks(np.arange(0.6, 0.73, 0.01))
    # if column_names[i] == 'mf_AUPR':
    #     ax.set_yticks(np.arange(0.27, 0.45, 0.01))
    # if column_names[i] == 'bp_Smin':
    #     ax.set_yticks(np.arange(51, 54, 0.2))
    # if column_names[i] == 'cc_Smin':
    #     ax.set_yticks(np.arange(11.8, 14, 0.2))
    # if column_names[i] == 'mf_Smin':
    #     ax.set_yticks(np.arange(15.5, 17.3, 0.2))








    # if 'Fmax' in column_names[i] or 'AUPR' in column_names[i]:
    #     ax.set_yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.01, 0.01))
    #     # ax.set_yticks(np.arange(min(y) - 0.01, max(y) + 0.01, 0.01))
    # elif 'Smin' in column_names[i]:
    #     ax.set_yticks(np.arange(np.floor(y_min), np.ceil(y_max) + 0.2, 0.2))
    #     # ax.set_yticks(np.arange(min(y) - 0.2, max(y) + 0.2, 0.2))


# 调整子图之间的间距和布局
plt.tight_layout()

# 调整图像分辨率
fig.set_dpi(300)

# 显示图像
plt.show()
print('lol')



'''
Single_plot::::

# 读取csv文件
data = pd.read_csv('CrossValidation_ColumnFilters.csv')

# 取第一列为横坐标
x = data.iloc[:, 0]
x = np.log2(x)

# 取包含"bp"的列，提取数据并作图
# bp_Fmax, bp_AUPR, cc_Fmax, cc_AUPR, mf_Fmax, mf_AUPR
for i in range(1, data.shape[1]):
    if "mf_Smin" in data.columns[i]:  # bp/cc/mf_Fmax/AUPR/Smin
        y = data.iloc[:, i]
        title = data.columns[i]
        marker = random.choice(markers)  # 随机选择一个标记
        plt.plot(x, y, label=title, linewidth=2, marker=marker, markersize=5)  # linestyle='-'

# 添加图例、横纵坐标标题
plt.legend(ncol=2, fontsize=8, loc='lower left')  # 调整图例字体大小 upper left, lower
plt.xlabel('log2(Filter)', fontsize=14)  # 调整横坐标标题字体大小
plt.xticks(range(3, 18, 1))  # range(3, 16, 1)
plt.ylabel('Value', fontsize=14)  # 调整纵坐标标题字体大小


# 调整图像分辨率和尺寸
fig = plt.gcf()
fig.set_size_inches(10, 8)  # 调整图像尺寸
fig.set_dpi(300)  # 调整图像分辨率

plt.show()
print('lol')
'''



'''
第一部分读取数据生成csv时候使用的chatGPT写代码输入prompt

写一个python代码，用于读取数据和生成csv表格：

root_path为第一层，其包含100个不同名称的文件夹，这100个文件夹为第二层，
文件名称包含Kernel和Filter:
DeepSS2GO_Kernel1024_Filter1024  
DeepSS2GO_Kernel128_Filter8     
DeepSS2GO_Kernel256_Filter64   
DeepSS2GO_Kernel4_Filter512     
DeepSS2GO_Kernel64_Filter4
DeepSS2GO_Kernel1024_Filter128   
kernel 取值范围 [4, 8, 16, 32, 64, 128, 256, 512, 1024] 
filter 取值范围 [4, 8, 16, 32, 64, 128, 256, 512, 1024] 

第二层的文件夹DeepSS2GO_Kernel64_Filter4包含的results文件夹为第三层
results文件夹中有多个txt文件为第四层，txt文件夹名称中包含[bp, cc, mf]

Python代码，读取每一个txt文件，提取其中的Fmax和AUPR信息，生成一个csv表格。

### 第一列是filter不同数值
csv表格第一列为第一层文件夹名称中filter后的数值，后面的每6列为第一层文件夹名中kernel后数值相同的一组，直到循环完所有的kernel数值。每一组6列的名称格式如下：
kernel4_bp_Fmax	kernel4_bp_AUPR	kernel4_cc_Fmax	kernel4_cc_AUPR	kernel4_mf_Fmax	kernel4_mf_AUPR

kernel8_bp_Fmax	kernel8_bp_AUPR	kernel8_cc_Fmax	kernel8_cc_AUPR	kernel8_mf_Fmax	kernel8_mf_AUPR

csv主体内容为对应的filter, kernel 和 Fmax，AUPR数值
比如：
DeepSS2GO_Kernel4_Filter512文件夹中results里的*bp.txt的Fmax数值，就要放在对应的filter=512的那一行，kernel4_bp_Fmax这一列
DeepSS2GO_Kernel16_Filter64文件夹中results里的*cc.txt的AUPR数值，就要放在对应的filter=64的那一行，kernel16_cc_AUPR这一列
DeepSS2GO_Kernel1024_Filter32文件夹中results里的*mf.txt的Fmax数值，就要放在对应的filter=32的那一行，kernel1024_mf_Fmax这一列

下面是一个txt文件的格式：
Length of test set: 662  threshold: 0.14
Smin: 53.814
Fmax: 0.353
AUPR: 0.278

'''



############################################################################
###### 第3部分，count_epochs，如果epochs很高，说明lr设置的太小 ###################
############################################################################



# 初始化一个字典用于记录每个filter和kernel对应的training.csv的行数
row_counts = {}
for f in filter_range:
    row_counts[f] = {}
    for k in kernel_range:
        row_counts[f][k] = 0

print(row_counts)

# 遍历第一层文件夹
for folder in os.listdir(root_path):
    folder_path = os.path.join(root_path, folder)
    if os.path.isdir(folder_path):
        # 如果是第二层文件夹，判断文件夹名是否符合要求
        if "Kernel" in folder and "Filter" in folder:
            # 获取filter和kernel的值
            filter_val = int(folder.split("_")[2][6:])
            kernel_val = int(folder.split("_")[1][6:])
            # 如果filter和kernel的值都在取值范围内，进一步遍历第三层文件夹
            if filter_val in filter_range and kernel_val in kernel_range:
                data_folder_path = os.path.join(folder_path, "data")
                if os.path.isdir(data_folder_path):
                    training_csv_path = os.path.join(data_folder_path, "training.csv")
                    if os.path.isfile(training_csv_path):
                        # 如果training.csv存在，计算行数
                        with open(training_csv_path, "r") as f:
                            row_count = sum(1 for line in f) - 1
                        row_counts[filter_val][kernel_val] = row_count
                        # print(training_csv_path)
                        # print(row_count)

# 生成csv表格
csv_path = "Count_Epochs.csv"  # 修改为你想输出的csv文件路径
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    header = ["Filter"] + [f"Kernel{k}" for k in kernel_range]
    writer.writerow(header)
    for f in filter_range:
        row = [f] + [row_counts[f][k] for k in kernel_range]
        writer.writerow(row)





'''
chatgpt promot

写一个python代码，用于读取数据和生成csv表格：

root_path为第一层，其包含100个不同名称的文件夹，这100个文件夹为第二层，
文件名称包含Kernel和Filter:
DeepSS2GO_Kernel1024_Filter1024  
DeepSS2GO_Kernel128_Filter8     
DeepSS2GO_Kernel256_Filter64   
DeepSS2GO_Kernel4_Filter512     
DeepSS2GO_Kernel64_Filter4
DeepSS2GO_Kernel1024_Filter128   
kernel 取值范围 [4, 8, 16, 32, 64, 128, 256, 512, 1024] 
filter 取值范围 [4, 8, 16, 32, 64, 128, 256, 512, 1024] 

第二层的文件夹DeepSS2GO_Kernel64_Filter4包含的data文件夹为第三层
data文件夹中有多个文件为第四层，其中一个是training.csv

Python代码，数一数training.csv文件有多少行（行数-1），生成一个csv表格。

第一列是filter不同数值
csv表格第一列为第一层文件夹名称中filter后的数值，后面每一列的title为kernel+数值，直到循环完所有的kernel数值
比如Kernel2, Kernel4, Kernel8, Kernel16等等

csv主体内容为对应的training.csv的（行数-1）
比如：
DeepSS2GO_Kernel4_Filter512文件夹中data/里的training.csv的行数，就要放在对应的filter=512的那一行，kernel4这一列
DeepSS2GO_Kernel16_Filter文件夹中data/里的training.csv的行数，就要放在对应的filter=的那一行，kernel16这一列
DeepSS2GO_Kernel1024_Filter32文件夹中data/里的training.csv的行数，就要放在对应的filter=32的那一行，kernel1024这一列

下面是一个training.csv文件的格式：此时的（行数-1）=12
epoch,loss,val_loss
0,0.1217421442270279,0.1006925031542778
1,0.1020296961069107,0.09912042319774628
2,0.09929966926574707,0.09760283678770065
3,0.09636150300502777,0.09594736248254776
4,0.09295417368412018,0.09538202732801437
5,0.08851592242717743,0.09393340349197388
6,0.08323285728693008,0.09353611618280411
7,0.07778183370828629,0.09308309853076935
8,0.07309440523386002,0.09356129914522171
9,0.0690111294388771,0.09367935359477997
10,0.0653148964047432,0.09353107959032059

输出csv文件名称为：epoch_counting.csv

'''



