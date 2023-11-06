# 这只是个放在“主文件夹”的例子，后面会在每个“子文件夹”生成独特的 params_local

params_local = {
    # from params_global_constant
    'device_ids': [0, 1],
    'run_step1_SplitTrainTest_Terms': 'T',
    'run_step2_Train': 'T',
    'run_step3_Test': 'T',
    'run_step4_pkl2fa': 'F',
    'run_step5_Diamond4CrossSpecies': 'F',
    'run_step6_FindAlpha': 'F',
    'run_step7_EvaluateAlpha': 'F',
    'run_step8_PredictAlpha': 'F',
    'aa_ss': 'ss8',
    'train_data': 'CAFA3',
    'test_data': 'CAFA3',
    'dir0': 'test_TrainECOLI_TestYEAST_aa_test/',
    'path_base': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/',
    'path_pub_data': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/',
    'go_file': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/pub_data/go.obo',
    'GOMinRepeat': 50,
    'TrainTestRatio': 0.95,
    'TrainValidRatio': 0.9,
    'batch_size': 32,
    'MAXLEN': 1024,
    'FC_depth': 0,
    'learning_rate': 0.0003,
    'epochs': 50,
    'load_pretrained_model': 0,
    'load_pretrained_model_addr': 'data/model_checkpoint.pth',
    'EarlyStopping_patience': 6,
    'PROT_LETTER_aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    'PROT_LETTER_ss8': ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B'],
    'PROT_LETTER_ss3': ['C', 'E', 'H'],

    # split: params_global_dynamic
    'kernels': [1],
    'filters': [1],
    'onts': 'all',
}


