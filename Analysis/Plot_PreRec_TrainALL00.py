import pandas as pd
from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
# from step0_TrainTestSetting_global import *

# README: 提取 results/PR_bp/cc/mf*pkl，画 Pression-recall曲线

# pre setting:
lw = 2

######### change following !!!!!!!!! ############

# path_base = '/home/fsong/work/py_proj/prot_algo/DeepSS2GO/'   # Keras 结果
path_base = '/home/fsong/work/py_proj/prot_algo/'  # Pytorch 结果




fpath_1 = path_base + 'DeepSS2GO/output/'

# fdir0aa = 'DeepGOGOGO_SPOT1DLM_TrainALL_TestALL_aa'  # ALL
# fdir0ss8 = 'DeepGOGOGO_SPOT1DLM_TrainALL_TestALL_ss8'  # ALL

fdir1 = 'test_TrainHUMAN_TestHUMAN_aa/DeepSS2GO_Kernel8_Filter65536_Ontsall'
fdir2 = 'test_TrainHUMAN_TestHUMAN_ss8/DeepSS2GO_Kernel48_Filter32768_Ontsall'

fdir3 = 'test_TrainHUMAN_TestARATH_aa/DeepSS2GO_Kernel16_Filter16384_Ontsall'
fdir4 = 'test_TrainHUMAN_TestARATH_ss8/DeepSS2GO_Kernel56_Filter65536_Ontsall'





# fdir5 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestMOUSE_aa'
# fdir6 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestMOUSE_ss8'
# fdir7 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestYEAST_aa'
# fdir8 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestYEAST_ss8'
# fdir9 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestSCHPO_aa'
# fdir10 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestSCHPO_ss8'
# fdir11 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestRAT_aa'
# fdir12 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestRAT_ss8'
# fdir13 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestECOLI_aa'
# fdir14 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestECOLI_ss8'
# fdir15 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestDROME_aa'
# fdir16 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestDROME_ss8'
# fdir17 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestCAEEL_aa'
# fdir18 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestCAEEL_ss8'
# fdir19 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestMYCTU_aa'
# fdir20 = 'DeepGOGOGO_SPOT1DLM_TrainHUMAN_TestMYCTU_ss8'





ont = 'mf'  # change among: bp, cc, mf  change !!!!
alpha = '1.00'  # change among: 0.00, 0.50, 1.00 change !!!!!

#####################################




# pre-setting
color = 'black'
linestyle = '--'


plt.figure()
def PR_plot(fpath_1, fdir, ont, alpha, color=color, linestyle=linestyle):
    df = pd.read_pickle('%s%s/results_alpha/PR_%s_%s.pkl' % (fpath_1, fdir, ont, alpha))
    x = df['recalls']
    y = df['precisions']
    aupr = np.trapz(y, x)
    plt.plot(x, y, color=color, linestyle=linestyle, lw=lw, label=f'{fdir[20:]} (AUPR = {aupr:0.2f})')  # fdir[20:] = 从20位往后截取名字


# plt change !!!!!!!!!!!!!!!
# PR_plot(fpath_1, fdir0aa, ont, alpha, 'black', '-')  # DeepGOPlus 复现
# PR_plot(fpath_1, fdir0ss8, ont, alpha, 'black', '--')  # DeepGOPlus 复现

PR_plot(fpath_1, fdir1, ont, alpha, 'black', '-')  # aa
PR_plot(fpath_1, fdir2, ont, alpha, 'black', '--')  # ss8
PR_plot(fpath_1, fdir3, ont, alpha, 'red', '-')  # red
PR_plot(fpath_1, fdir4, ont, alpha, 'red', '--')  # red
# PR_plot(fpath_1, fdir5, ont, alpha, 'orange', '-')
# PR_plot(fpath_1, fdir6, ont, alpha, 'orange', '--')
# PR_plot(fpath_1, fdir7, ont, alpha, 'yellow', '-')
# PR_plot(fpath_1, fdir8, ont, alpha, 'yellow', '--')
# PR_plot(fpath_1, fdir9, ont, alpha, 'green', '-')
# PR_plot(fpath_1, fdir10, ont, alpha, 'green', '--')

# PR_plot(fpath_1, fdir11, ont, alpha, 'blue', '-')
# PR_plot(fpath_1, fdir12, ont, alpha, 'blue', '--')
# PR_plot(fpath_1, fdir13, ont, alpha, 'purple', '-')
# PR_plot(fpath_1, fdir14, ont, alpha, 'purple', '--')
# PR_plot(fpath_1, fdir15, ont, alpha, 'grey', '-')
# PR_plot(fpath_1, fdir16, ont, alpha, 'grey', '--')
# PR_plot(fpath_1, fdir17, ont, alpha, 'pink', '-')
# PR_plot(fpath_1, fdir18, ont, alpha, 'pink', '--')
# PR_plot(fpath_1, fdir19, ont, alpha, 'magenta', '-')
# PR_plot(fpath_1, fdir20, ont, alpha, 'magenta', '--')



plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title('(' + ont + ') ')
plt.legend(loc="upper right", prop={'size': 6})


plt.show()


