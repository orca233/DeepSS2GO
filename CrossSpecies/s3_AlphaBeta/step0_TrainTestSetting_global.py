import os
# import configparser
# import numpy as np

# step0_TrainTestSetting_global.py 中包含 params_global_dynamic & params_global_constant 两个字典，
# 读取 params_global_dynamic，拆分list成元素，更新到 params_global_constant，
# 保存在 step0_TrainTestSetting_local.py 中的 params_local 字典

print('##################### step0 starting ##############################')

#################### path_base ################
## fpath  判断root根目录，看是HPC还是lab的workstation，设置path_base=绝对路径
dir_sustech_hpc = '/scem/work/songfu/py_proj/prot_algo/DeepSS2GO/'  # change -----------------
lab_linux3090 = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'  # change -----------------

path_base = ''
# 判断该路径是否存在，决定path_base是采用HPC集群路径，或是实验室的Liao_lab路径
if os.path.exists(dir_sustech_hpc):
    path_base = dir_sustech_hpc
elif os.path.exists(lab_linux3090):
    path_base = lab_linux3090
else:
    print('HOLY SHIT, NO --path_base-- available !!!')


# # 相对路径path_base
# # 获取当前文件夹路径
# current_directory = os.getcwd()
# # 选择上推两层的文件夹路径
# path_base = os.path.abspath(os.path.join(current_directory, '..', '..'))
# path_base = path_base + '/'
# print('path_base = ', path_base)


################################################
############# params_global_constant ###########
################################################

params_global_constant = {
    # 硬件参数
    'device_ids': [0, 1],  # 单一gpu运算：str='cuda:3',   多GPU并行: list=[0, 1, 2, 3]， CPU计算: str='cpu'

    ##### 下面这几行内容，是打算分开计算 Train & Test
    # 是否运行哪几个step，
    'run_step1_SplitTrainTest_Terms': 'T',
    'run_step2_Train': 'T',  # 如果 train = T，则 不加载 pretrained model，需要预训练
    'run_step3_Test': 'T',  #
    'run_step4_pkl2fa': 'F',  #
    'run_step5_Diamond4CrossSpecies': 'F',  #
    'run_step6_FindAlpha': 'F',
    'run_step7_EvaluateAlpha': 'F',  # 评价三大指标，Fmax, AUPR, Smin，分别用带alpha的和不带alpha的
    'run_step8_PredictAlpha': 'F',

    # 下面这几个可能会和 global_dynamic的互换：
    'aa_ss': 'aa',  # aa, ss8, ss3 三个选项  ['aa', 'ss8']
    'train_data': 'ECOLI',  #  只写具体物种 HUMAN, MOUSE, ARATH。。。，不用写ALL00， ['HUMAN', 'MOUSE']
    'test_data': 'YEAST',  # ['HUMAN', 'MOUSE']

    # 文件夹
    # dir0 = 是“非变量，不参与循环”       dir1 = 变量，参与循环
    'dir0': 'test_TrainECOLI_TestYEAST_aa_test/',  # output/dir0/ 同一批次实验的root文件夹， change -----

    'path_base': path_base,
    'path_pub_data': path_base + 'pub_data/',
    'go_file': path_base + 'pub_data/go.obo',
    # 'path_redundancy': path_base + 'redundancy/',

    # model 软件训练参数
    'GOMinRepeat': 50,  # 最少GO term出现次数, GO Min Repeat - gominre
    'TrainTestRatio': 0.95,  # (train + valid) / test  - traintestratio
    'TrainValidRatio': 0.90,  # train / valid
    'batch_size': 32,
    'MAXLEN': 1024,  # 最长的氨基酸数目，DeepGOPlus为2000，这里因为SPOT1DLM最长1021，所以设为: 1024
    'FC_depth': 0,  # 替换为所需的全连接层深度, original = 3 -- dense_depth
    'learning_rate': 3e-4,
    'epochs': 50,  # 记得改回 10000 ！！！！！！！！！！！!!!!!!!!!!!!!!! 否则kernel=4, filter=4的时候，epoch跑到100也跑不完

    'load_pretrained_model': 0,  # 是否加载 已经训练好的模型
    'load_pretrained_model_addr': 'data/model_checkpoint.pth',  # 要修改

    'EarlyStopping_patience': 6,
    'PROT_LETTER_aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    'PROT_LETTER_ss8': ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B'],  # no E in aa_letter
    'PROT_LETTER_ss3': ['C', 'E', 'H'],  # C = C+S+T    E = E+B     H = H+G+I



}


################################################
############# params_global_dynamic ###########
################################################

# 把list拆开分散到每个子文件夹，kernels & filters，可作为loop，也可整体作为List，kernel目前最大也就1024了，因为受到esm限制
# 这里必须写list，否则在run_step1_2.py中就会循环HUMAN的每一个字母了
# split_list = ['train_data', 'test_data']  # 要拆分哪些list

params_global_dynamic = {
    # ont 一定参与循环：可以是 ['all'] -- 训练所有三种terms，或 ['bp', 'cc', 'mf'] -- 逐一训练
    # 在 param_local.py 只有一个str: all/bp/cc/mf 四选一
    'onts': ['all'],  # onts 有几种形式： ['bp'], ['cc'], ['mf'], ['gominre_trxte'], ['gominre], ['bp', 'cc', 'mf']  # 旧版本有 ['all'], ['gominre']只针对>50条件1

    'kernels': [8, 16],
    'filters': [16, 32]   # 不需要变成 tuple

    # 顺序
    # kernels = [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 256]  # 512
    # filters = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]  # 131072

    # 倒序 reverse:
    # 'kernels': [256, 128, 120, 112, 104, 96, 88, 80, 72, 64, 56, 48, 40, 32, 24, 16, 8]  # 512
    # 'filters': [65536, 32768, 16384, 8192, 4096, 2048, 1024, 512, 256, 128, 64, 32, 16]   # 131072

    # 并行
    # 'kernels': [[8, 8], [16, 16], [24, 24], [131072, 131072]],
    # 'filters': [[32, 32], [64, 64]]   # 不需要变成 tuple

}

print('step0 done ')


### 一些尝试：###

###### TrainHUMAN_TestHUMAN_KernelX_FilterY  ######
# 单核下：cuda:3  VS  [3]，区别不大
# kernel=8, filter=262,144,   HPC V100，双核不可算。            lab_linux_3090，单核不可算，双核不可算，四核不可算。即使onts只算bp/cc/mf也无用

# kernel=128, filter=16,384,                                       lab_linux_3090，单核(epoch=4.5min)，双核(epoch=2min)，四核(epoch=2min)
# kernel=128, filter=32,768,                                       lab_linux_3090，单核(epoch=8.5min)，双核(epoch=4min)，四核(epoch=4min)
# kernel=128, filter=65,536,                                       lab_linux_3090，单核不可算，双核(epoch=7min)，四核(epoch=7min)
# kernel=128, filter=131,072,   HPC_V100，单核不可算，双核可算?       lab_linux_3090，单核不可算，双核不可算，四核(epoch=17min) ，等同于***

# kernel=256, filter=16,384,                                       lab_linux_3090，单核(epoch=7min)，双核(epoch=3min)，四核(epoch=3min)
# kernel=256, filter=32,768,                                       lab_linux_3090，单核(epoch=14min)，双核(epoch=5.5min)，四核(epoch=5.5min)
# kernel=256, filter=65,536,                                       lab_linux_3090，单核不可算，双核(epoch=12min)，四核(epoch=12min)
# kernel=256, filter=131,072,   HPC V100，                         lab_linux_3090，单核不可算，双核不可算，四核不可算

# kernel=[128, 128], filter=[65536, 65536]                       lab_linux_3090，单核不可算，双核不可算，四核(epoch=14min)，等同于***
# kernel=[128, 128], filter=[131072, 131072]                      lab_linux_3090，单核不可算，双核不可算，四核不可算
# kernel=[64, 64, 64, 64], filter=[65536, 65536, 65536, 65536]    lab_linux_3090，单核不可算，双核不可算，四核不可算

# kernel=[8, 8], filter=[131072, 131072]             HPC_V100，34核，epoch=1h20min       lab_linux_3090，cpu可算，36核，epoch=2h,
# kernel=[16, 16], filter=[131072, 131072]，         HPC_V100，34核，epoch=1h40min
# kernel=[24, 24], filter=[131072, 131072]，         HPC_V100，34核，epoch=2h10min

###### TrainALL00_TestALL00_KernelX_FilterY  ######
# kernel=128, filter=32,768,    HPC_V100，单核(epoch=45min)，双核(epoch=27min)     lab_linux_3090，单核(epoch=min)，双核(epoch=min)，四核(epoch=min)
# kernel=128, filter=65,536,    HPC_V100，单核(epoch=1.5h)，双核(epoch=53min)      lab_linux_3090，单核不可算，双核(epoch=min)，四核(epoch=min)
# kernel=128, filter=131,072,   HPC_V100，单核不可算，双核不可算 ！！！                lab_linux_3090，单核不可算，双核不可算，四核(epoch=min)

# kernel=256, filter=32,768,    HPC_V100，单核(epoch=1h12min)，双核(epoch=40min)      lab_linux_3090，单核(epoch=min)，双核(epoch=min)，四核(epoch=min)
# kernel=256, filter=65,536,    HPC_V100，单核(epoch=2.5h)，双核(epoch=1h20min)       lab_linux_3090，单核不可算，双核(epoch=min)，四核(epoch=min)
# kernel=256, filter=131,072,   HPC V100，单核不可算，双核不可算                     lab_linux_3090，单核不可算，双核不可算，四核不可算?

# HPC GPU 直接计算 K=8, F=131072。  可以计算  K=[[4, 4], [8, 8], [12, 12], [16, 16], [24, 24], [32, 32], [40, 40], [48, 48]], F=[[65536, 65536]]， K56以上崩了, epoch=1h


###### TrainHUMAN_TestHUMAN_DeepGOPlus ######
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter16384_16384_...*8，lab_linux_3090，单核不可算，双核不可算，四核不可算
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter8192_8192_...*8，lab_linux_3090，单核不可算，双核不可算，四核(epoch=12min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter4096_4096_...*8，lab_linux_3090，单核不可算，双核(epoch=8min)，四核(epoch=6min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter2048_2048_...*8，lab_linux_3090，单核不可算，双核(epoch=4min)，四核(epoch=3min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter1024_1024_...*8，lab_linux_3090，单核不可算，双核(epoch=2min)，四核(epoch=Xmin)


###### TrainALL00_TestALL00_DeepGOPlus ######
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter16384_16384_...*8，lab_linux_3090，单核不可算，双核不可算，四核不可算
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter8192_8192_...*8，lab_linux_3090，单核不可算，双核不可算，四核不可算
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter4096_4096_...*8，lab_linux_3090，单核不可算，双核(epoch=47min)，四核(epoch=48min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter2048_2048_...*8，lab_linux_3090，单核不可算，双核(epoch=24min)，四核(epoch=26min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter1024_1024_...*8，lab_linux_3090，单核不可算，双核(epoch=13min)，四核(epoch=13min)
# Kernel8_16_24_32_40_48_56_64_72_80_88_96_104_112_120_128_Filter512_512_...*8，lab_linux_3090，单核不可算，双核(epoch=Xmin)，四核(epoch=7min)


'''
kernels = [8, 16, 24], filters = [32, 64]，
经过run循环后，生成 6个 dir1：K8F32, K8F64, K16F32, K16F64, K24F32, K24F64

kernels = [[8, 8], [16, 16], [24, 24]], filters': [[32, 32], [64, 64]]
经过run循环后，生成 6个 dir1: K[8, 8]F[32, 32],   K[8, 8]F[64, 64],   K[16, 16]F[32, 32], ......
对于K[8, 8]F[32, 32]，带入模型变成：K8--F64
对于K[8, 8]F[32768, 32768]，带入模型变成：K8--F65536 ！！！！！！
'''


# 判断是否运行TrainALL00_TestALL00，或者TrainX_TestY, =0不运行，=1运行。？？？？？？？？？这俩有啥区别来着？？？？？？？？？？？？
# 答案：如果运行ALL00，要考虑的是，不要TrainALL00_TestHUMAN  ！！ 在 两个list里会交叉选，
# 比如 train_data = [ALL00, HUMAN, MOUSE], test_data = [ALL00, ECOLI]，就会出现 TrainALL00_TestECOLI






