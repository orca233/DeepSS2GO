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


def process_folder(folder_path):
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


def main():
    base_path = "/data1/fsong/py_proj/prot_algo/DeepSS2GO/output/test_CrossSpecies/"
    target_folder = "test_TrainHUMAN_TestARATH_aa"
    folder_path = os.path.join(base_path, target_folder)

    if os.path.exists(folder_path):
        results = process_folder(folder_path)

        # Creating the DataFrame with specific row indices
        data = {
            'metrics': ['bp_Fmax', 'bp_AUPR', 'bp_Smin', 'cc_Fmax', 'cc_AUPR', 'cc_Smin', 'mf_Fmax', 'mf_AUPR',
                        'mf_Smin'],
            'score': [results['bp']['Fmax'], results['bp']['AUPR'], results['bp']['Smin'],
                      results['cc']['Fmax'], results['cc']['AUPR'], results['cc']['Smin'],
                      results['mf']['Fmax'], results['mf']['AUPR'], results['mf']['Smin']],
            'kernel_filter': [results['bp']['kernel_fmax'], results['bp']['kernel_aupr'], results['bp']['kernel_smin'],
                              results['cc']['kernel_fmax'], results['cc']['kernel_aupr'], results['cc']['kernel_smin'],
                              results['mf']['kernel_fmax'], results['mf']['kernel_aupr'], results['mf']['kernel_smin']]
        }

        df = pd.DataFrame(data).transpose()
        df.columns = df.iloc[0]  # Set first row as column header
        df = df.drop(df.index[0])  # Drop the initial row as it's now the header

        print(df)
        print('lol\n')
        print(df['bp_AUPR'])
        print(df.loc['score'])

    else:
        print(f"Folder not found: {folder_path}")


if __name__ == "__main__":
    main()
