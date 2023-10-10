import os
params_global_constant = {
    # 文件夹
    # 下面这几个可能会和 global_dynamic的互换：
    'aa_ss': ['aa'],  # aa, ss8, ss3 三个选项  ['aa', 'ss8']
    'onts': ['bp'],  # ['bp', 'cc', 'mf']
    'train_data': ['HUMAN', 'MOUSE'],  #  只写具体物种 HUMAN, MOUSE, ARATH。。。，不用写ALL00， ['HUMAN', 'MOUSE']
    'test_data': ['HUMAN', 'ARATH']  # ['HUMAN', 'MOUSE']
}

split_list = ['train_data', 'test_data']

params_global_dynamic = {
    'train_data': ['HUMAN', 'MOUSE'],
    'test_data': ['HUMAN', 'ARATH'],
    'aa_ss': ['aa'],
    'kernels': [16, 24],
    'filters': [128, 512],
    'onts': ['bp']
}

for kernel_temp in params_global_dynamic['kernels']:   # ----- change -----
    for filter_temp in params_global_dynamic['filters']:   # ----- change -----

        step0_TrainTestSetting_local_fpath = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO_Pytorch/output/test_TrainHUMAN_TestHUMAN_KernelXFilterY_aassZ_OntK/test_local.py'

        with open(step0_TrainTestSetting_local_fpath, "w") as local_file:
            local_file.write("params_local = {\n")

            # 不变的数据 global_constant
            for key, value in params_global_constant.items():
                local_file.write(
                    f"    '{key}': {repr(value)},\n")  # original: (f"    '{key}': {repr(value)},\n")，会把原始格式保留存到local.py里, e.g. [[16, 24]]




            # local_file.write(f"    'kernels': [{value}],\n")


            m = 0
            for key, value in params_global_dynamic.items():
                print(m)
                m += 1



                if key in split_list:  # 在split_list列表，拆！  从list转成单个元素
                    local_file.write(f"    '{key}': {value},\n")
                else:  # 在split_list列表，不拆！
                    local_file.write(f"    '{key}': {repr(value)},\n")  # original: (f"    '{key}': {repr(value)},\n")

            local_file.write("}\n")

