'''
README:
用于


小tip:
include_kernel_filter = True  # Set to False to exclude kernel_filter
这里可决定是否在 DataFrame 中包含 kernel_filter 数据，还是仅包含 score数据

'''

import os
import re
import pandas as pd

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

def process_folder(folder_path, include_kernel_filter):
    results = {
        'bp': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''},
        'cc': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''},
        'mf': {'Fmax': 0, 'AUPR': 0, 'Smin': float('inf'), 'kernel_fmax': '', 'kernel_aupr': '', 'kernel_smin': ''}
    }

    for subfolder in os.listdir(folder_path):
        if re.match(r'DeepSS2GO_Kernel(\d+)_Filter(\d+)_Ontsall', subfolder):
            kernel, filter = re.findall(r'Kernel(\d+)_Filter(\d+)', subfolder)[0]
            results_path = os.path.join(folder_path, subfolder, "results_alpha")
            for file in os.listdir(results_path):
                if file.endswith('.txt'):
                    category = file.split('_')[3]  # category = bp/cc/mf
                    if category not in results:
                        continue

                    metrics = extract_metrics(os.path.join(results_path, file))
                    kernel_filter = [int(kernel), int(filter)]
                    if 'Fmax' in metrics and metrics['Fmax'] > results[category]['Fmax']:
                        results[category]['Fmax'] = metrics['Fmax']
                        results[category]['kernel_fmax'] = kernel_filter
                    if 'AUPR' in metrics and metrics['AUPR'] > results[category]['AUPR']:
                        results[category]['AUPR'] = metrics['AUPR']
                        results[category]['kernel_aupr'] = kernel_filter
                    if 'Smin' in metrics and metrics['Smin'] < results[category]['Smin']:
                        results[category]['Smin'] = metrics['Smin']
                        results[category]['kernel_smin'] = kernel_filter

    return results

def create_dataframe(data, include_kernel_filter, train, test, folder):
    df = pd.DataFrame(data).transpose()
    df.columns = df.iloc[0]
    df = df.drop(df.index[0])
    df['folder'] = folder
    df['Train'] = train
    df['Test'] = test

    if not include_kernel_filter:
        df = df.drop(columns=['kernel_filter'])

    return df

def main():
    base_path = "/data1/fsong/py_proj/prot_algo/DeepSS2GO/output/test_CrossSpecies/"
    all_dfs_aa = []
    all_dfs_ss8 = []
    include_kernel_filter = False  # Set to True / False to exclude kernel_filter  这里可决定是否在 DataFrame 中包含 kernel_filter 数据，还是仅包含 score数据

    folder_k = 0
    for folder in os.listdir(base_path):
        match = re.match(r'test_Train(\w+)_Test(\w+)_(aa|ss8)', folder)
        if match:
            print('folder_k = ', folder_k)
            folder_k += 1
            train, test, suffix = match.groups()
            folder_path = os.path.join(base_path, folder)
            if os.path.exists(folder_path):
                results = process_folder(folder_path, include_kernel_filter)
                data = {
                    'metrics': ['bp_Fmax', 'bp_AUPR', 'bp_Smin', 'cc_Fmax', 'cc_AUPR', 'cc_Smin', 'mf_Fmax', 'mf_AUPR', 'mf_Smin'],
                    'score': [results['bp']['Fmax'], results['bp']['AUPR'], results['bp']['Smin'],
                              results['cc']['Fmax'], results['cc']['AUPR'], results['cc']['Smin'],
                              results['mf']['Fmax'], results['mf']['AUPR'], results['mf']['Smin']],
                    'kernel_filter': [results['bp']['kernel_fmax'], results['bp']['kernel_aupr'], results['bp']['kernel_smin'],
                                      results['cc']['kernel_fmax'], results['cc']['kernel_aupr'], results['cc']['kernel_smin'],
                                      results['mf']['kernel_fmax'], results['mf']['kernel_aupr'], results['mf']['kernel_smin']]
                }

                df = create_dataframe(data, include_kernel_filter, train, test, folder)

                if suffix == 'aa':
                    all_dfs_aa.append(df)
                elif suffix == 'ss8':
                    all_dfs_ss8.append(df)

    # Concatenate and sort all dataframes
    final_df_aa = sort_and_save(all_dfs_aa, 'CrossSpecies_results_aa.csv')
    final_df_ss8 = sort_and_save(all_dfs_ss8, 'CrossSpecies_results_ss8.csv')
    print("AA Results:\n", final_df_aa)
    print("SS8 Results:\n", final_df_ss8)

def sort_and_save(dfs, filename):
    final_df = pd.concat(dfs)
    order = ['ARATH', 'ECOLI', 'HUMAN', 'MOUSE', 'MYCTU', 'YEAST']
    final_df['Train'] = pd.Categorical(final_df['Train'], categories=order, ordered=True)
    final_df['Test'] = pd.Categorical(final_df['Test'], categories=order, ordered=True)
    final_df = final_df.sort_values(by=['Train', 'Test'])
    final_df.to_csv(filename, index=False)
    return final_df

if __name__ == "__main__":
    main()
