import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the datasets
# path_CrossSpecies = '/mnt/data/'  # /mnt/data/
path_CrossSpecies = '/data1/fsong/py_proj/prot_algo/DeepSS2GO/Analysis/CrossSpecies_HeatMap/'  # /mnt/data/
df_aa = pd.read_csv(path_CrossSpecies + 'CrossSpecies_results_aa.csv')
df_ss8 = pd.read_csv(path_CrossSpecies + 'CrossSpecies_results_ss8.csv')

# Species for train and test
species = ["ARATH", "ECOLI", "HUMAN", "MOUSE", "MYCTU", "YEAST"]

# Metrics to choose from
metrics = ['Fmax', 'AUPR', 'Smin']

# Ask user for metric choice
# metric_choice = input("Choose the metric for the heatmap ('Fmax', 'AUPR', 'Smin'): ")  # 提供一个选项
metric_choice = 'Fmax'  # 'Fmax', 'AUPR', 'Smin' 三选一

# Check if the choice is valid
if metric_choice not in metrics:
    raise ValueError(f"Invalid metric choice. Please choose from {metrics}.")

# Define a function to create a heatmap data matrix
def create_heatmap_matrix(df, feature):
    matrix = []
    for test_species in species:
        row = []
        for train_species in species:
            value = df.loc[
                (df['Train'] == train_species) & (df['Test'] == test_species),
                feature
            ].values[0]
            row.append(value)
        matrix.append(row)
    return matrix

# Create heatmap matrices for the chosen metric
heatmap_data = {
    f'bp_{metric_choice}_aa': create_heatmap_matrix(df_aa, f'bp_{metric_choice}'),
    f'bp_{metric_choice}_ss8': create_heatmap_matrix(df_ss8, f'bp_{metric_choice}'),
    f'cc_{metric_choice}_aa': create_heatmap_matrix(df_aa, f'cc_{metric_choice}'),
    f'cc_{metric_choice}_ss8': create_heatmap_matrix(df_ss8, f'cc_{metric_choice}'),
    f'mf_{metric_choice}_aa': create_heatmap_matrix(df_aa, f'mf_{metric_choice}'),
    f'mf_{metric_choice}_ss8': create_heatmap_matrix(df_ss8, f'mf_{metric_choice}')
}

# Plot the heatmaps
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle(metric_choice, fontsize=16)  # 图片最上方大标题

# Define the colormap for each set of heatmaps
cmap_dict = {
    f'bp_{metric_choice}_aa': 'Greens',
    f'bp_{metric_choice}_ss8': 'Greens',
    f'cc_{metric_choice}_aa': 'Blues',
    f'cc_{metric_choice}_ss8': 'Blues',
    f'mf_{metric_choice}_aa': 'Oranges',
    f'mf_{metric_choice}_ss8': 'Oranges'
}

# Titles for the heatmaps
titles = [f'bp_{metric_choice}_aa', f'cc_{metric_choice}_aa', f'mf_{metric_choice}_aa',
          f'bp_{metric_choice}_ss8', f'cc_{metric_choice}_ss8', f'mf_{metric_choice}_ss8']

# Increase font size for the annotations in the heatmaps
sns.set(font_scale=1.2)

# Reorder the axes for plotting as per the new layout
axes_order = [0, 2, 4, 1, 3, 5]
axes = [axes.flat[i] for i in axes_order]

for ax, title in zip(axes, titles):
    sns.heatmap(heatmap_data[title], ax=ax, cmap=cmap_dict[title], annot=True, fmt=".3f",
                xticklabels=species, yticklabels=species)
    ax.set_title(title.replace('_', ' '))
    ax.set_xlabel('Train')
    ax.set_ylabel('Test')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for main title
plt.savefig(path_CrossSpecies + f'CrossSpecies_HeatMap_{metric_choice}.png')
plt.show()

'''
1）如果metric_choice是Fmax或AUPR，则图片中的数字为3位小数。如果是Smin，则图片中数字为2位小数。
2）将bp的两个图aa和ss8上下排列在最左边，cc两个图的aa和ss8上下排列在中间，mf两个图的aa和ss8上下排列在最右边
'''
