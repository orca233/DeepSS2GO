import os


### fpath
dir_SUSTech_HPC = '/scem/work/songfu/py_proj/prot_algo/DeepSS2GO/'  # change ----------------- original: DeepSS2GO_Pytorch
dir_liao_lab = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'  # change ----------------- original: DeepSS2GO_Pytorch

if os.path.exists(dir_SUSTech_HPC):
    path_base = dir_SUSTech_HPC
elif os.path.exists(dir_liao_lab):
    path_base = dir_liao_lab
else:
    print('NO path_base avail')

# # 获取当前文件夹路径
# current_directory = os.getcwd()
# # 选择上推两层的文件夹路径
# path_base = os.path.abspath(os.path.join(current_directory, '..', '..'))
# path_base = path_base + '/'
# print('path_base = ', path_base)


# path_Prot_T5_XL_UniRef50 = '/scem/work/songfu/py_proj/prot_data/Prot_T5_XL_UniRef50'  # HPC change ---------------
path_Prot_T5_XL_UniRef50 = '/home/fsong/work/py_proj/prot_data/Prot_T5_XL_UniRef50'  # local lab_linux
# path_esm1b = '/home/fsong/work/py_proj/prot_data/esm'  # change ---------------

path_pub_data = path_base + 'pub_data/'
path_SPOT1DLM_checkpoints = path_base + 'pub_data/SPOT1DLM_checkpoints/'
path_redundancy = path_base + 'redundancy/'



# SPOT1DLM 生成琐碎文件redundancy位置：
# /scem/work/songfu/prot_algo/DeepSS2GO/redundancy/sSPOT1DLM_inputs
# save_path_npy = path_redundancy + 'SPOT1DLM_inputs/'  # for step3 & 4, esm/prottrans生成*npy
# os.system('mkdir -p %s' % save_path_npy)
save_path_npy = path_base + 'redundancy/SPOT1DLM_inputs_new/'  # for step3 & 4, esm/prottrans生成*npy
os.system('mkdir -p %s' % save_path_npy)

# /scem/work/songfu/prot_algo/DeepSS2GO/redundancy/sSPOT1DLM_results
# save_path_csv = path_redundancy + 'SPOT1DLM_results/'  # for step5 inference生成的 *csv文件
# os.system('mkdir -p %s' % save_path_csv)
save_path_csv = path_redundancy + 'SPOT1DLM_results_new/'  # for step5 inference生成的 *csv文件
os.system('mkdir -p %s' % save_path_csv)

print('step0 done')

