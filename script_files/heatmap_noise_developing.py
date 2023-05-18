import os
import math
import pandas as pd
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib import rcParams
import pysam
import numpy as np
import matplotlib.colors as mcolors
import glob

def get_average_depth(bam_file, reference_name):
    alignment_file = pysam.AlignmentFile(bam_file, "rb")
    total_depth = 0
    reference_length = alignment_file.get_reference_length(reference_name)
    for column in alignment_file.pileup(reference_name):
        total_depth += column.n
    average_depth = total_depth / reference_length
    return average_depth

def shannon_entropy(base_counts):
    entropy = 0
    total_count = sum(base_counts.values())

    for base, count in base_counts.items():
        probability = count / total_count
        entropy -= probability * math.log2(probability)

    return entropy

def noise_scores(bam_file, reference_name):
    alignment_file = pysam.AlignmentFile(bam_file, "rb")
    reference_length = alignment_file.get_reference_length(reference_name)
    noise_scores = [0] * reference_length

    for column in alignment_file.pileup(reference_name):
        base_counts = Counter()
        for read in column.pileups:
            if not read.is_del and not read.is_refskip:
                base = read.alignment.query_sequence[read.query_position]
                base_counts[base] += 1

        noise_scores[column.reference_pos] = shannon_entropy(base_counts)

    return noise_scores

def save_noise_scores_to_csv(sample, bam_file, reference_name, output_csv):
    scores = noise_scores(bam_file, reference_name)
    data = {
        "sample": [sample] * len(scores),
        "reference": [reference_name] * len(scores),
        'position': list(range(len(scores))),
        'noise': scores
    }
    df = pd.DataFrame(data)
    return df

bam_dir = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/bam"
bam_files = glob.glob(f"{bam_dir}/*.bam")
samples = [os.path.basename(bam_file).split('.')[0] for bam_file in bam_files]

reference_name = "A_HA_H1"  # Reference name
output_csv = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/"

dfs = []
for sample, bam_file in zip(samples, bam_files):
    df = save_noise_scores_to_csv(sample, bam_file, reference_name, output_csv)
    dfs.append(df)

print(dfs)

# Define the color palette
'''
colors = ['#a9d6e5', '#468faf', '#2a6f97', '#01497c', '#012a4a']
bins = [-np.inf, 10, 100, 200, np.inf]
cmap = sns.color_palette(colors, n_colors=len(bins) - 1, as_cmap=True)
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(bins, cmap.N, clip=True)
'''

# Choose the 'viridis' colormap
cmap = plt.get_cmap('viridis')

# Concatenate the DataFrames into one combined DataFrame
combined_df = pd.concat(dfs)

# Pivot the combined DataFrame to have 'sample' as index, 'position' as columns, and 'noise' as values
pivoted_df = combined_df.pivot(index='sample', columns='position', values='noise')

# Create the heatmap using the style provided
fig, ax = plt.subplots(figsize=(10, 8))
sns.heatmap(pivoted_df, cmap=cmap, annot=False, linewidths=0, ax=ax)


# Customize the chart appearance
ax.set_title('H1N1 Gene Depth', fontsize=18, color='#255F85')
ax.set_xlabel('')
fig.patch.set_facecolor('#F5F5DC')
ax.set_facecolor('#F5F5DC')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(axis='x', colors='#255F85', labelsize=6)
ax.tick_params(axis='y', colors='#255F85', labelsize=10)
border = Rectangle((0, 0), 1, 1, transform=ax.transAxes, fill=False, edgecolor='#255F85', linewidth=2.5)
ax.add_patch(border)

# Show the chart
plt.show()