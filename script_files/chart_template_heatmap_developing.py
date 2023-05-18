import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib import rcParams
import numpy as np
import matplotlib.colors as mcolors

# Data
data = pd.read_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results/new_influenza_pipline_final.csv')
df = pd.DataFrame(data)
df = df.dropna(subset=['Subtype'])

df_H1N1 = df[df['Subtype'].str.contains('H1N1')]
df_H3N2 = df[df['Subtype'].str.contains('H3N2')]
df_VIC = df[df['Subtype'].str.contains('VIC')]
columns_to_keep_H1N1 = ['sample', 'Avg_Depth_A_HA_H1', 'Avg_Depth_A_H1_MP', 'Avg_Depth_A_NA_N1', 'Avg_Depth_A_H1_NS', 'Avg_Depth_A_H1_NP', 'Avg_Depth_A_H1_PA', 'Avg_Depth_A_H1_PB1', 'Avg_Depth_A_H1_PB2']
columns_to_keep_H3N2 = ['sample', 'Avg_Depth_A_HA_H3', 'Avg_Depth_A_H3_MP', 'Avg_Depth_A_NA_N2', 'Avg_Depth_A_H3_NS', 'Avg_Depth_A_H3_NP', 'Avg_Depth_A_H3_PA', 'Avg_Depth_A_H3_PB1', 'Avg_Depth_A_H3_PB2']
columns_to_keep_VIC = ['sample', 'Avg_Depth_B_VIC_HA', 'Avg_Depth_B_VIC_NA']

df_H1N1 = df_H1N1[columns_to_keep_H1N1]
df_H3N2 = df_H3N2[columns_to_keep_H3N2]
df_VIC = df_VIC[columns_to_keep_VIC]

# Fill NaN values with 0
df_H1N1.fillna(0, inplace=True)
df_H3N2.fillna(0, inplace=True)
df_VIC.fillna(0, inplace=True)

# Convert all columns (except 'sample') to float data type and round them to integer
cols_to_convert = [col for col in df_H1N1.columns if col != 'sample']
df_H1N1[cols_to_convert] = df_H1N1[cols_to_convert].astype(float).round().astype(int)
cols_to_convert = [col for col in df_H3N2.columns if col != 'sample']
df_H3N2[cols_to_convert] = df_H3N2[cols_to_convert].astype(float).round().astype(int)
cols_to_convert = [col for col in df_VIC.columns if col != 'sample']
df_VIC[cols_to_convert] = df_VIC[cols_to_convert].astype(float).round().astype(int)

# Create a list of the dataframes and column names for the three charts
charts_data = [(df_H1N1[columns_to_keep_H1N1], 'H1N1 Gene Depth'),
               (df_H3N2[columns_to_keep_H3N2], 'H3N2 Gene Depth'),
               (df_VIC[columns_to_keep_VIC], 'VIC Gene Depth')]

# Remove "Avg_Depth_" from the column names on the y-axis
df_H1N1.columns = df_H1N1.columns.str.replace('Avg_Depth_', '')
df_H3N2.columns = df_H3N2.columns.str.replace('Avg_Depth_', '')
df_VIC.columns = df_VIC.columns.str.replace('Avg_Depth_', '')



#HEATMAP TEMPLATE
# Set the Bodoni 72 font (assuming it's installed on your system)
rcParams['font.family'] = 'Bodoni 72'
rcParams['font.weight'] = 'bold'

# Define the color palette
colors = ['#a9d6e5', '#468faf', '#2a6f97', '#01497c', '#012a4a']
bins = [-np.inf,10,100,200, np.inf]
cmap = sns.color_palette(colors, n_colors=len(bins)-1, as_cmap=True)
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bins, cmap.N, clip=True)


# Set up the Seaborn style to resemble Vox charts
sns.set_style("whitegrid", {'grid.linestyle': '--'})
sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})


# Create a heatmap
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(df_H1N1.set_index('sample').T, cmap=cmap, norm=norm, annot=False, fmt="d", linewidths=0.7, cbar=False, ax=ax, square=True)

# Set the x-axis label to an empty string
ax.set_xlabel('')

# Customize chart appearance
ax.set_title('H1N1 Gene Depth', fontsize=18, color='#255F85')

# Set background color to light beige
fig.patch.set_facecolor('#F5F5DC')
ax.set_facecolor('#F5F5DC')

# Remove default spines
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Change tick colors
ax.tick_params(axis='x', colors='#255F85', labelsize=6)
ax.tick_params(axis='y', colors='#255F85', labelsize=10)

# Add a thicker border around the chart area
border = Rectangle((0, 0), 1, 1, transform=ax.transAxes,fill=False, edgecolor='#255F85', linewidth=2.5)
ax.add_patch(border)

# Show the chart
plt.show()
