import os
# import configparser
# import numpy as np


print('##################### step0 starting ##############################')


dir_sustech_hpc = '/scem/work/songfu/py_proj/prot_algo/DeepSS2GO/'
lab_linux3090 = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'

path_base = ''
if os.path.exists(dir_sustech_hpc):
    path_base = dir_sustech_hpc
elif os.path.exists(lab_linux3090):
    path_base = lab_linux3090
else:
    print('SH, NO --path_base-- available !!!')


params_global_constant = {

    'device_ids': [0, 1],
    'aa_ss': 'aa',  # aa, ss8, ss3   ['aa', 'ss8']
    'train_data': 'MOUSE',
    'test_data': 'HUMAN',
    'dir0': 'test_TrainX_TestY_aa/',

    'run_step1_SplitTrainTest_Terms': 'T',
    'run_step2_Train': 'T',
    'run_step3_Test': 'T',
    'run_step4_pkl2fa': 'T',

    'run_step7.1_EvaluateWithoutAlpha': 'T',
    'run_step8.1_PredictWithoutAlpha': 'F',

    'run_step5_Diamond4CrossSpecies': 'F',
    'run_step6_FindAlpha': 'F',
    'run_step7_EvaluateAlpha': 'F',
    'run_step8_PredictAlpha': 'F',

    'path_base': path_base,
    'path_pub_data': path_base + 'pub_data/',
    'go_file': 'data/go.obo',


    # model
    'GOMinRepeat': 50,
    'TrainTestRatio': 0.95,
    'TrainValidRatio': 0.90,
    'batch_size': 32,
    'MAXLEN': 1024,
    'FC_depth': 0,
    'learning_rate': 3e-4,
    'epochs': 50,

    'load_pretrained_model': 0,
    'load_pretrained_model_addr': 'data/model_checkpoint.pth',

    'EarlyStopping_patience': 6,
    'PROT_LETTER_aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    'PROT_LETTER_ss8': ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B'],
    'PROT_LETTER_ss3': ['C', 'E', 'H'],  # C = C+S+T    E = E+B     H = H+G+I

}


############# params_global_dynamic ###########

params_global_dynamic = {
    'onts': ['all'],
    'kernels': [8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128],
    'filters': [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]


}

print('step0 done ')

