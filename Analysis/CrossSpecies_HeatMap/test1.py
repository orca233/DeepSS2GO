import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Load the datasets
df_aa = pd.read_csv('/mnt/data/CrossSpecies_results_aa.csv')
df_ss8 = pd.read_csv('/mnt/data/CrossSpecies_results_ss8.csv')

# Species for train and test
species = ["ARATH", "ECOLI", "HUMAN", "MOUSE", "MYCTU", "YEAST"]

# Define a function to create a heatmap data matrix
def create_heatmap_matrix_corrected(df, feature):
    matrix = []
    for test_species in species:
        row = []
        for train_species in species:
            # Extract the value for the specific train-test pair
            value = df.loc[
                (df['Train'] == train_species) & (df['Test'] == test_species),
                feature
            ].values[0]
            row.append(value)
        matrix.append(row)
    return matrix

# Create heatmap matrices for the corrected features
heatmap_data_corrected = {
    'bp_Fmax_aa': create_heatmap_matrix_corrected(df_aa, 'bp_Fmax'),
    'bp_Fmax_ss8': create_heatmap_matrix_corrected(df_ss8, 'bp_Fmax'),
    'cc_Fmax_aa': create_heatmap_matrix_corrected(df_aa, 'cc_Fmax'),
    'cc_Fmax_ss8': create_heatmap_matrix_corrected(df_ss8, 'cc_Fmax'),
    'mf_Fmax_aa': create_heatmap_matrix_corrected(df_aa, 'mf_Fmax'),
    'mf_Fmax_ss8': create_heatmap_matrix_corrected(df_ss8, 'mf_Fmax')
}

# Plot the corrected heatmaps
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Define the colormap for each set of heatmaps
cmap_dict_corrected = {
    'bp_Fmax_aa': 'Greens',
    'bp_Fmax_ss8': 'Greens',
    'cc_Fmax_aa': 'Blues',
    'cc_Fmax_ss8': 'Blues',
    'mf_Fmax_aa': 'Oranges',
    'mf_Fmax_ss8': 'Oranges'
}

# Titles for the heatmaps
titles = ['bp_Fmax_aa', 'cc_Fmax_aa', 'mf_Fmax_aa', 'bp_Fmax_ss8', 'cc_Fmax_ss8', 'mf_Fmax_ss8']

# Increase font size for the annotations in the heatmaps
sns.set(font_scale=1.2)  # Adjust this to increase or decrease font size

# Reorder the axes for plotting as per the new layout
axes_order = [0, 2, 4, 1, 3, 5]
axes = [axes.flat[i] for i in axes_order]

for ax, title in zip(axes, titles):
    sns.heatmap(heatmap_data_corrected[title], ax=ax, cmap=cmap_dict_corrected[title], annot=True, fmt=".2f",
                xticklabels=species, yticklabels=species)
    ax.set_title(title.replace('_', ' ').title())
    ax.set_xlabel('Train')
    ax.set_ylabel('Test')

plt.tight_layout()
plt.savefig('/mnt/data/heatmap_figures_corrected.png')
plt.show()
