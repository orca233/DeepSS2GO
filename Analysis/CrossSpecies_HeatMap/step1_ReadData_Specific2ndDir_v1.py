import os
import re


def extract_metrics(file_path):
    metrics = {}
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split(": ")
            if len(parts) == 2:
                key, value = parts
                try:
                    metrics[key] = float(value)
                except ValueError:
                    continue  # Skip lines where conversion to float fails
    return metrics


def process_folder(folder_path):
    results = {
        'bp': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''},
        'cc': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''},
        'mf': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''}
    }

    for subfolder in os.listdir(folder_path):
        # folder_path=/data1/fsong/py_proj/prot_algo/DeepSS2GO/output/test_CrossSpecies/test_TrainHUMAN_TestARATH_aa/
        if re.match(r'DeepSS2GO_Kernel\d+_Filter\d+_Ontsall', subfolder):
            # print('22222 subfolder = ', subfolder)
            results_path = os.path.join(folder_path, subfolder, "results_alpha")
            for file in os.listdir(results_path):
                if file.endswith('.txt'):
                    category = file.split('_')[3]  # category = bp/cc/mf
                    # print('333 category=', category)
                    if category not in results:
                        continue  # Skip if category is not in the results dictionary

                    metrics = extract_metrics(os.path.join(results_path, file))
                    if 'Fmax' in metrics and metrics['Fmax'] > results[category]['Fmax']:
                        results[category]['Fmax'] = metrics['Fmax']
                        results[category]['kernel_fmax'] = subfolder
                    if 'AUPR' in metrics and metrics['AUPR'] > results[category]['AUPR']:
                        results[category]['AUPR'] = metrics['AUPR']
                        results[category]['kernel_aupr'] = subfolder
                    if 'Smin' in metrics and metrics['Smin'] < results[category]['Smin']:
                        results[category]['Smin'] = metrics['Smin']
                        results[category]['kernel_smin'] = subfolder

    return results



def main():
    base_path = "/data1/fsong/py_proj/prot_algo/DeepSS2GO/output/test_CrossSpecies/"
    target_folder = "test_TrainHUMAN_TestARATH_aa"  # Change this to your specific folder
    folder_path = os.path.join(base_path, target_folder)

    if os.path.exists(folder_path):
        print('111111111folder_path = ', folder_path)
        results = process_folder(folder_path)
        for category in ['bp', 'cc', 'mf']:
            print(f"{category} Fmax: {results[category]['Fmax']} (Kernel/Filter: {results[category]['kernel_fmax']})")
            print(f"{category} AUPR: {results[category]['AUPR']} (Kernel/Filter: {results[category]['kernel_aupr']})")
            print(f"{category} Smin: {results[category]['Smin']} (Kernel/Filter: {results[category]['kernel_smin']})")
    else:
        print(f"Folder not found: {folder_path}")

if __name__ == "__main__":
    main()





'''
在linux服务器上有一组数据。
一级文件夹地址为：/data1/fsong/py_proj/prot_algo/DeepSS2GO/output/test_CrossSpecies/，
二级文件夹为：test_TrainM_TestN_aa/，或 test_TrainM_TestN_ss8/
M和N是[ARATH, ECOLI, HUMAN, MOUSE, MYCTU, YEAST]中的一个。
注意：有些二级文件夹的命名不规则，则忽略。

三级文件夹为：DeepSS2GO_KernelX_FilterY_Ontsall/
X和Y是两个数字，
X取值范围：[8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128]
Y取值范围：[16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]

四级文件夹为：results_alpha/
四级文件夹内有3个txt文件，分别为
Fmax_AUPR_Smin_bp_1.00.txt
Fmax_AUPR_Smin_cc_1.00.txt
Fmax_AUPR_Smin_mf_1.00.txt

在每一个txt中，格式如下：
Smin: 21.772
Fmax: 0.304
AUPR: 0.232


写一个python程序，遍历某一个二级文件夹（比如test_TrainHUMAN_TestARATH_aa/）下的所有三级文件夹内的四级文件夹，读取所有txt的内容。
输出9个数值，[bp_Fmax，cc_Fmax, mf_Fmax, bp_AUPR, cc_AUPR, mf_AUPR, bp_Smin, cc_Smin, mf_Smin]

bp_Fmax为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_bp_1.00.txt内Fmax最大的数字。
cc_Fmax为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_cc_1.00.txt内Fmax最大的数字。
mf_Fmax为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_mf_1.00.txt内Fmax最大的数字。


bp_AUPR为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_bp_1.00.txt内AUPR最大的数字。
cc_AUPR为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_cc_1.00.txt内AUPR最大的数字。
mf_AUPR为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_mf_1.00.txt内AUPR最大的数字。

bp_Smin为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_bp_1.00.txt内Smin最小的数字。
cc_Smin为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_cc_1.00.txt内Smin最小的数字。
mf_Smin为三级文件夹不同Kernel和Filter下的四级文件夹内，所有Fmax_AUPR_Smin_mf_1.00.txt内Smin最小的数字。

除了这9个数值，还要给出每个数值取得最大(Fmax/AUPR)或最小值(Smin)时，所对应的Kernel和Filter的X/Y值



最终输出结果为dataframe格式。
第一行是9个名称，第二行是数值value，第三行是对应的[kernel, filter]

df['TrainHUMAN_TestARATH'] = \
[
['bp_Fmax', 'bp_AUPR', 'bp_Smin', 'cc_Fmax', 'cc_AUPR', 'cc_Smin', 'mf_Fmax', 'mf_AUPR', 'mf_Smin'],
[0.339, 0.265, 24.59, 0.745, 0.769, 5.825, 0.43, 0.379, 6.611],
[[16, 65536], [16, 65536], [16, 65536], [24, 65536], [24, 65536], [24, 65536], [16, 16384], [8, 8192], [16, 16384]]
]
'''





'''
输出格式为：
bp Fmax: 0.339 (Kernel/Filter: DeepSS2GO_Kernel16_Filter65536_Ontsall)
bp AUPR: 0.265 (Kernel/Filter: DeepSS2GO_Kernel16_Filter65536_Ontsall)
bp Smin: 24.59 (Kernel/Filter: DeepSS2GO_Kernel16_Filter65536_Ontsall)
cc Fmax: 0.745 (Kernel/Filter: DeepSS2GO_Kernel24_Filter65536_Ontsall)
cc AUPR: 0.769 (Kernel/Filter: DeepSS2GO_Kernel24_Filter65536_Ontsall)
cc Smin: 5.825 (Kernel/Filter: DeepSS2GO_Kernel24_Filter65536_Ontsall)
mf Fmax: 0.43 (Kernel/Filter: DeepSS2GO_Kernel16_Filter16384_Ontsall)
mf AUPR: 0.379 (Kernel/Filter: DeepSS2GO_Kernel8_Filter8192_Ontsall)
mf Smin: 6.611 (Kernel/Filter: DeepSS2GO_Kernel16_Filter16384_Ontsall)

修改代码，使得输出为dataframe格式如下：
df['TrainHUMAN_TestARATH'] = \
[
['bp_Fmax', 'bp_AUPR', 'bp_Smin', 'cc_Fmax', 'cc_AUPR', 'cc_Smin', 'mf_Fmax', 'mf_AUPR', 'mf_Smin'],
[0.339, 0.265, 24.59, 0.745, 0.769, 5.825, 0.43, 0.379, 6.611],
[[16, 65536], [16, 65536], [16, 65536], [24, 65536], [24, 65536], [24, 65536], [16, 16384], [8, 8192], [16, 16384]]
]


'''
