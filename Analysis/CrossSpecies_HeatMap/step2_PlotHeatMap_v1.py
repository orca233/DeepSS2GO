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
    'A': create_heatmap_matrix_corrected(df_aa, 'bp_Fmax'),
    'B': create_heatmap_matrix_corrected(df_ss8, 'bp_Fmax'),
    'C': create_heatmap_matrix_corrected(df_aa, 'cc_Fmax'),
    'D': create_heatmap_matrix_corrected(df_ss8, 'cc_Fmax'),
    'E': create_heatmap_matrix_corrected(df_aa, 'mf_Fmax'),
    'F': create_heatmap_matrix_corrected(df_ss8, 'mf_Fmax')
}

# Plot the corrected heatmaps
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Define the colormap for each set of heatmaps
cmap_dict_corrected = {'A': 'Greens', 'B': 'Greens', 'C': 'Blues', 'D': 'Blues', 'E': 'Oranges', 'F': 'Oranges'}

for ax, (label, data_matrix) in zip(axes.flat, heatmap_data_corrected.items()):
    sns.heatmap(data_matrix, ax=ax, cmap=cmap_dict_corrected[label], annot=True, fmt=".2f",
                xticklabels=species, yticklabels=species)
    ax.set_title(f"Heatmap of Data = {label.lower()}_Fmax")
    ax.set_xlabel('X-axis (Train)')
    ax.set_ylabel('Y-axis (Test)')

plt.tight_layout()
plt.savefig('/mnt/data/heatmap_figures_corrected.png')
plt.show()
