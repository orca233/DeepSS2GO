# README： 交叉验证是在做对于下面两组参数kernels & filters 一对一的交叉炼丹，看相互影响
# 这里面的train_data & test_data设置一个值，就不来回切换了

# 运行命令： 后台 nohup python run_CrossValidation0_diffKernelsFilters_step1step2.py &

import os
import sys
from step0_TrainTestSetting_global import *  # 这个是从 step0_global.py中提取信息，生成 step0_local.py
# 循环读取每个配置项并添加到字典中
# import configparser
from typing import Union

###### 在DeepSS2GO/中新建 output/dir0/dir1/(data/ or results/...)
current_dir = os.getcwd()  # 后面要用

os.chdir(path_base)  # 退回到上一层，DeepSS2GO/
os.system('mkdir -p output')
os.chdir('output')  # 进入DeepSS2GO/output/

# dir0同一批测试,只有文件夹名字，不含绝对路径path，是一系列测试的local_root, e.g. dir0 = test_diffKernelsFilters/
os.system('mkdir -p %s' % params_global_constant['dir0'])  # 在 output/里新建的dir0/
dir0_path = '%soutput/%s' % (path_base, params_global_constant['dir0'])  # 是指向dir0的绝对路径
# os.system('mkdir -p %s_redundancy' % dir0_path)  # 平行在output/建立一个dir0_path_redundancy/文件夹，放model.h5等超大的文件


####################################################################
######################### TrainX_TestY ########################
####################################################################

# 循环所有，根据实际情况修改：这里的list可以只有一个元素，也可以有很多元素（会逐一生成对应文件夹）

# for train_data_temp in params_global_dynamic['train_data']:
#     for test_data_temp in params_global_dynamic['test_data']:
#         for aa_ss_temp in params_global_dynamic['aa_ss']:
#             for kernel_temp in params_global_dynamic['kernels']:
#                 for filter_temp in params_global_dynamic['filters']:
#                     for ont_temp in params_global_dynamic['onts']:
#                         # CreateDirections_CreateParamslocal_RunSteps()


### 方法一： 手动拆分 某些list，比如 kernel & filter
# 这个循环根据实际情况进行修改
for kernel_temp in params_global_dynamic['kernels']:   # ----- change -----
    for filter_temp in params_global_dynamic['filters']:   # ----- change -----
        for ont_temp in params_global_dynamic['onts']:  # ----- change -----

            # 第一步： 生成一系列文件夹: output/dir0/dir1/      dir1里包含 data/ & results/ 等 ###########
            # 判断 kernel_temp 是一个 list=[16, 16, 16] 还是 int = 16。  如果是前者，则输出16_16_16
            def process_input(x):
                if isinstance(x, list):
                    y = '_'.join(str(num) for num in x)
                else:
                    y = x
                return y

            # print('111111')
            # print(kernel_temp)
            # print(filter_temp)

            dir1 = f'DeepSS2GO_Kernel{process_input(kernel_temp)}_Filter{process_input(filter_temp)}_Onts{ont_temp}/'  # ----- change -----
            print('dir1 = ', dir1)
            dir1_path = os.path.join(dir0_path, dir1)
            os.makedirs(dir1_path, exist_ok=True)  # exist_ok 类似于 mkdir -p，已存在则跳过。创建 output/dir0/dir1/文件夹
            # os.system('cp %s*.py %s' % (path_base + 's2_TrainTest/', dir1_path))  # cp s2_TrainTest/*py 到 dir1/

            # print('111111111111')
            # print(current_dir)
            os.system(f'cp {current_dir}/*py {dir1_path}')  # cp 当前文件夹s2_TrainTest_xxx/*py 到 dir1/
            # os.system(f'cp -rf {current_dir}/data {dir1_path}')  ################### 仅针对于CAFA3！！！！！！！！！！！！！！！#########################  更新的CAFA3是不需要的

            os.chdir(dir1_path)  # 进入该文件夹 output/dir0/dir1
            os.system('mkdir -p data')  # dir1/下新建data/，有了也无妨
            os.system('mkdir -p results')
            # os.chdir(dir0_path)  # 退回到output/文件夹


            # 第二步：生成：dir1/step0_TrainTestSetting_local.py， 包含params_local字典
            # 包含params_global_constant的所有信息，
            # 并将params_global_dynamic 的list分裂后，把单个元素添加到params_local字典中
            step0_TrainTestSetting_local_fpath = os.path.join(dir1_path, "step0_TrainTestSetting_local.py")
            with open(step0_TrainTestSetting_local_fpath, "w") as local_file:
                local_file.write("params_local = {\n")
                local_file.write("    # from params_global_constant\n")

                # 不变的数据 global_constant
                for key, value in params_global_constant.items():
                    local_file.write(
                        f"    '{key}': {repr(value)},\n")  # original: (f"    '{key}': {repr(value)},\n")，会把原始格式保留存到local.py里, e.g. [[16, 24]]



                # "可能"拆分的数据，global_dynamic
                local_file.write("\n    # split: params_global_dynamic\n")

                if isinstance(kernel_temp, int):
                    local_file.write(f"    'kernels': [{kernel_temp}],\n")  # 这样输出的是list，而不是一个int
                else:
                    local_file.write(f"    'kernels': {kernel_temp},\n")  # 这样输出的是 单一list，而不是多重list，比如[[64, 64]]

                if isinstance(filter_temp, int):
                    local_file.write(f"    'filters': [{filter_temp}],\n")  # 这样输出的是list，而不是一个int
                else:
                    local_file.write(f"    'filters': {filter_temp},\n")  # 这样输出的是 单一list，而不是多重list，比如[[131072, 131072]]

                local_file.write(f"    'onts': {repr(ont_temp)},\n")  # 加上repr，这样输出的是一个str


                # for key, value in params_global_dynamic.items():
                    # if key in split_list:  # 在split_list列表，拆！  从list转成单个元素
                    #     # local_file.write(f"    '{key}': [{value}],\n")
                    #     local_file.write(f"    '{key}': {value},\n")
                    # else:  # 在split_list列表，不拆！
                    #     local_file.write(f"    '{key}': {repr(value)},\n")  # original: (f"    '{key}': {repr(value)},\n")

                local_file.write("}\n")


            # # 第三步：判断是否运行step 1 2 3
            os.chdir(dir1_path)  # 要在local运算
            if params_global_constant['run_step1_splittraintest_terms'] == 'T':  # 只运行step1的话，可以用来在output/生成一系列dir1文件夹，然后做具体的step2参数炼丹实验。
                os.system('python step1_SplitTrainTest_Terms_Pytorch.py')
            if params_global_constant['run_step2_TrainTest'] == 'T':
                os.system('python step2_TrainTest_Pytorch.py')
            if params_global_constant['run_step3_evaluate'] == 'T':
                os.system('python step3_Evaluate_Pytorch.py')  # 第三步分运行，否则顺延太耽误时间

            os.chdir(dir0_path)  # 回到dir0



print('run_TrainTest done')


















# # 方法二，split_list中的元素为将要从params_global_dynamic进行分裂的key。不知道一共有几个key。。。这是个很麻烦的事情
# # 比如下面：只分裂train_data & test_data 参与循环，但是kernels & filters保持为list迭代到 locol.py中
#
# split_list = ['train_data', 'test_data']
#
# params_global_dynamic = {
#     'train_data': ['HUMAN', 'MOUSE'],  #  只写具体物种 HUMAN, MOUSE, ARATH。。。，不用写ALL00， ['HUMAN', 'MOUSE']
#     'test_data': ['HUMAN', 'ARATH'],  # ['HUMAN', 'MOUSE']
#     'aa_ss': ['aa'],  # aa, ss8, ss3 三个选项  ['aa', 'ss8']
#     'kernels': [16, 24],  # [4, 8, 16, 24, 40, 48, 56, 72, 80, 88, 96, 104, 112, 120]
#     'filters': [128, 512],  # [16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768]  # 不需要变成 tuple
#     'onts': ['bp']  # ['bp', 'cc', 'mf']
# }

# #### 下面这是通过递归方法，自动把所有的在split_list的key进行循环，
#
# split_list = ['train_data', 'test_data']
#
# params_global_dynamic = {
#     'train_data': ['HUMAN', 'MOUSE'],
#     'test_data': ['HUMAN', 'ARATH'],
#     'aa_ss': ['aa'],
#     'kernels': [16, 24],
#     'filters': [128, 512],
#     'onts': ['bp']
# }
#
# def recursive_loop(split_list, params_global_dynamic, current_index, values):
#     if current_index == len(split_list):
#         # 达到循环的最后一层，执行你需要的操作
#         print(values)
#         return
#
#     current_key = split_list[current_index]
#     current_values = params_global_dynamic[current_key]
#
#     for value in current_values:
#         updated_values = values + [value]
#         print('1111111')
#         print('updated_values=', updated_values)
#         recursive_loop(split_list, params_global_dynamic, current_index + 1, updated_values)
#
# # 从第一个索引开始递归循环
# recursive_loop(split_list, params_global_dynamic, 0, [])

'''
这段代码定义了一个名为 recursive_loop 的递归函数。
它接受 split_list、params_global_dynamic、loop_index 和 loop_values 作为参数。
loop_index 表示当前循环嵌套的层级，loop_values 是一个列表，用于存储每个循环层级的值。
在每次递归调用中，我们检查当前的 loop_index 是否达到了 split_list 的长度。
如果是，则说明已经完成了所有循环层级，我们可以处理循环的结果。
否则，我们遍历 params_global_dynamic[split_list[loop_index]] 的值，并将其添加到 loop_values 中。
然后，我们递归调用 recursive_loop，将 loop_index 增加 1，并继续下一层的循环。
完成一次递归后，我们将之前添加的值从 loop_values 中弹出，以便进行下一次循环。
通过这种递归的方式，可以根据 split_list 的长度动态确定循环的嵌套次数，并获取每个循环层级的值。
你可以根据需要修改 recursive_loop 函数中的代码，以执行你需要的操作。
'''

