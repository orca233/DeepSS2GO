import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the datasets
path_CrossSpecies = '/data1/fsong/py_proj/prot_algo/DeepSS2GO/Analysis/CrossSpecies_HeatMap/'  # /mnt/data/
# path_CrossSpecies = '/mnt/data/'
df_aa = pd.read_csv(path_CrossSpecies + 'CrossSpecies_results_aa.csv')
df_ss8 = pd.read_csv(path_CrossSpecies + 'CrossSpecies_results_ss8.csv')

# Species for train and test
species = ["ARATH", "ECOLI", "HUMAN", "MOUSE", "MYCTU", "YEAST"]

# Metrics to choose from
metrics = ['Fmax', 'AUPR', 'Smin']

# Set metric choice
metric_choice = 'AUPR'  # 'Fmax', 'AUPR', 'Smin' 三选一 ！！！！！！！！！！！！！！！

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
fig.suptitle(metric_choice, fontsize=16)

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

# Determine the number format based on the metric choice
if metric_choice in ['Fmax', 'AUPR']:
    fmt = ".3f"  # 3 decimal places for Fmax and AUPR
elif metric_choice == 'Smin':
    fmt = ".2f"  # 2 decimal places for Smin

# Increase font size for the annotations in the heatmaps
sns.set(font_scale=1.2)

# Create a list to order the axes properly (3 columns by 2 rows)
axes_order = [(i, j) for i in range(2) for j in range(3)]

# Reorder the axes for plotting as per the new layout
for ((row, col), title) in zip(axes_order, titles):
    sns.heatmap(heatmap_data[title], ax=axes[row, col], cmap=cmap_dict[title], annot=True, fmt=fmt,
                xticklabels=species, yticklabels=species)
    axes[row, col].set_title(title.replace('_', ' '))
    axes[row, col].set_xlabel('Train')
    axes[row, col].set_ylabel('Test')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for main title
heatmap_output_path = path_CrossSpecies + f'CrossSpecies_HeatMap_{metric_choice}.png'
plt.savefig(heatmap_output_path)
plt.show()

# # Return the path to the saved heatmap for the user to download
# heatmap_output_path
