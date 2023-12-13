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

# Set metric choice to 'Fmax' since we are using that data
metric_choice = 'Smin'

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

# Find min and max for each pair of heatmaps to set the color bar range

# original:
# bp_min = min(map(min, heatmap_data['bp_Fmax_aa'] + heatmap_data['bp_Fmax_ss8']))
# bp_max = max(map(max, heatmap_data['bp_Fmax_aa'] + heatmap_data['bp_Fmax_ss8']))
#
# cc_min = min(map(min, heatmap_data['cc_Fmax_aa'] + heatmap_data['cc_Fmax_ss8']))
# cc_max = max(map(max, heatmap_data['cc_Fmax_aa'] + heatmap_data['cc_Fmax_ss8']))
#
# mf_min = min(map(min, heatmap_data['mf_Fmax_aa'] + heatmap_data['mf_Fmax_ss8']))
# mf_max = max(map(max, heatmap_data['mf_Fmax_aa'] + heatmap_data['mf_Fmax_ss8']))


# FS change
bp_min = min(map(min, heatmap_data[f'bp_{metric_choice}_aa'] + heatmap_data[f'bp_{metric_choice}_ss8']))
bp_max = max(map(max, heatmap_data[f'bp_{metric_choice}_aa'] + heatmap_data[f'bp_{metric_choice}_ss8']))

cc_min = min(map(min, heatmap_data[f'cc_{metric_choice}_aa'] + heatmap_data[f'cc_{metric_choice}_ss8']))
cc_max = max(map(max, heatmap_data[f'cc_{metric_choice}_aa'] + heatmap_data[f'cc_{metric_choice}_ss8']))

mf_min = min(map(min, heatmap_data[f'mf_{metric_choice}_aa'] + heatmap_data[f'mf_{metric_choice}_ss8']))
mf_max = max(map(max, heatmap_data[f'mf_{metric_choice}_aa'] + heatmap_data[f'mf_{metric_choice}_ss8']))





# Set up the matplotlib figure and axes
fig, axes = plt.subplots(2, 3, figsize=(18, 12))
fig.suptitle(metric_choice, fontsize=16)

# Titles for the heatmaps
titles = [f'bp_{metric_choice}_aa', f'cc_{metric_choice}_aa', f'mf_{metric_choice}_aa',
          f'bp_{metric_choice}_ss8', f'cc_{metric_choice}_ss8', f'mf_{metric_choice}_ss8']

# Define the colormap for each set of heatmaps
cmap_dict = {
    'bp': 'Greens',
    'cc': 'Blues',
    'mf': 'Oranges'
}

# Determine the number format based on the metric choice
fmt = ".3f"
sns.set(font_scale=1.2)

# Create a list to order the axes properly (3 columns by 2 rows)
axes_order = [(i, j) for i in range(2) for j in range(3)]

# Plot each heatmap
for ((row, col), title) in zip(axes_order, titles):
    ax = axes[row, col]
    category = title.split('_')[0]
    vmin, vmax = eval(f'{category}_min'), eval(f'{category}_max')
    sns.heatmap(heatmap_data[title], ax=ax, cmap=cmap_dict[category], annot=True, fmt=fmt,
                xticklabels=species, yticklabels=species, vmin=vmin, vmax=vmax)
    ax.set_title(title.replace('_', ' '))
    ax.xaxis.set_ticks_position('top')  # Set the x-axis labels to the top
    ax.xaxis.set_label_position('top')  # Set the x-axis labels' position to the top
    ax.set_xlabel('Train')
    ax.set_ylabel('Test')

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust for main title
plt.show()
