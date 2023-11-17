# 这只是个放在“主文件夹”的例子，后面会在每个“子文件夹”生成独特的 params_local


params_local = {
    # from params_global_constant
    'device_ids': [0, 1],
    'run_step1_splittraintest_terms': 'T',
    'run_step2_TrainTest': 'T',
    'run_step3_evaluate': 'F',
    'dir0': 'test_TrainMYCTU_TestHUMAN_aa/',
    # 'path_base': '../../',  # 从output往前推两格
    'path_base': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/',

    # 'path_pub_data': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/pub_data/',
    # 'path_redundancy': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/redundancy/',
    # 'path_s2_TrainTest': '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/s2_TrainTest/',
    # 'go_file': '../../pub_data/go.obo',
    'GOMinRepeat': 50,
    'TrainTestRatio': 0.95,
    'TrainValidRatio': 0.9,
    'batch_size': 32,
    'MAXLEN': 1024,
    'FC_depth': 0,
    'learning_rate': 0.0003,
    'epochs': 50,
    'load_pretrained_model': 0,
    'load_pretrained_model_addr': 'xxx/model_checkpoint.pth',
    'EarlyStopping_patience': 6,
    'PROT_LETTER_aa': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    'PROT_LETTER_ss8': ['C', 'S', 'T', 'H', 'G', 'I', 'E', 'B'],
    'PROT_LETTER_ss3': ['C', 'E', 'H'],
    'aa_ss': 'aa',
    'train_data': 'ALL00',
    'test_data': 'ALL00',

    # split: params_global_dynamic
    'kernels': [16],
    'filters': [65536],
    'onts': 'all',

    # step5_FindAlpha
    'alpha_range': [25, 80, 1]

}


