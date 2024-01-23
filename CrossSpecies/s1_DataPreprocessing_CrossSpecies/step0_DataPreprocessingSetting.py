import os

### fpath
dir_SUSTech_HPC = '/scem/work/songfu/py_proj/prot_algo/DeepSS2GO/'  # change -----------------
dir_liao_lab = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'  # change -----------------

path_base = 'see following'
if os.path.exists(dir_SUSTech_HPC):
    path_base = dir_SUSTech_HPC
elif os.path.exists(dir_liao_lab):
    path_base = dir_liao_lab
else:
    print('NO path_base avail')


path_Prot_T5_XL_UniRef50 = '/scem/work/songfu/py_proj/prot_data/Prot_T5_XL_UniRef50'  # HPC change

path_pub_data = path_base + 'pub_data/'
path_SPOT1DLM_checkpoints = path_base + 'pub_data/SPOT1DLM_checkpoints/'
path_redundancy = path_base + 'redundancy/'

save_path_npy = path_redundancy + 'SPOT1DLM_inputs/'
os.system('mkdir -p %s' % save_path_npy)

save_path_csv = path_redundancy + 'SPOT1DLM_results/'
os.system('mkdir -p %s' % save_path_csv)

print('step0 done')

