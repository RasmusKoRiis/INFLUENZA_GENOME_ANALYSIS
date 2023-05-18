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

    # Add 0 for missing bases
    all_bases = {'A', 'C', 'G', 'T'}
    missing_bases = all_bases - set(base_counts.keys())
    for base in missing_bases:
        base_counts[base] = 0

    # Calculate the base proportions
    base_proportions = {base: count / total_count for base, count in base_counts.items()}

    # Sort the bases by their proportions in descending order
    sorted_bases = sorted(base_proportions.keys(), key=lambda base: base_proportions[base], reverse=True)

    # Check if the three bases that are not in the majority of reads have roughly the same proportion
    proportion_threshold = 0.1
    non_majority_proportions = [base_proportions[base] for base in sorted_bases[1:]]
    real_noise = all(abs(p1 - p2) <= proportion_threshold for p1 in non_majority_proportions for p2 in non_majority_proportions)

    # Calculate the Shannon entropy
    for base, count in base_counts.items():
        probability = count / total_count
        if probability > 0:  # Skip calculation for bases with probability 0
            entropy -= probability * math.log2(probability)

    return entropy, real_noise, base_proportions



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

def save_noise_scores_to_csv(sample, bam_file, reference_name, output_csv, fig, ax, min_base_proportion=0.1):
    scores, real_noise_flags, base_ratios = zip(*noise_scores(bam_file, reference_name))
    data = {
        "sample": [sample] * len(scores),
        "reference": [reference_name] * len(scores),
        'position': list(range(len(scores))),
        'noise': scores,
        'real_noise': real_noise_flags,
        'base_ratios': base_ratios
    }
    df = pd.DataFrame(data)

    df.to_csv(output_csv, index=False)

    bars = ax.bar(df['position'], df['noise'], color='#255F85', edgecolor='black', linewidth=1)

    threshold = 0.5
    for idx, bar in enumerate(bars):
        bar_height = bar.get_height()
        if bar_height > threshold:
            ax.scatter(bar.get_x() + bar.get_width() / 2, bar_height, color='red', marker='o', s=30, zorder=2)

    print(df)

def customize_and_display_charts(fig, axs, samples):
    rcParams['font.family'] = 'Bodoni 72'
    rcParams['font.weight'] = 'bold'

    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})

    for idx, ax in enumerate(axs.flat):
        if idx < len(samples):
            ax.set_title(samples[idx], fontsize=18, color='#255F85', pad=20)
        ax.set_xlabel('Position', fontsize=14, color='#255F85')
        ax.set_ylabel('Noise', fontsize=14, color='#255F85')
        ax.set_ylim(0, 1.0)

        fig.patch.set_facecolor('#F5F5DC')
        ax.set_facecolor('#F5F5DC')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        ax.tick_params(axis='x', colors='#255F85')
        ax.tick_params(axis='y', colors='#255F85')

        ax.yaxis.grid(False)

        border = Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                        fill=False, edgecolor='#255F85', linewidth=2.5)
        ax.add_patch(border)
    fig.savefig('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/combined_figure.png', dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())

if __name__ == "__main__":
    bam_directory = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/bam"
    reference_name = "A_HA_H1"
    output_csv = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/"

    bam_files = [f for f in os.listdir(bam_directory) if f.endswith('.bam')]
    samples = [os.path.splitext(os.path.basename(bam_file))[0] for bam_file in bam_files]

    n_bam_files = len(bam_files)
    nrows = math.ceil(n_bam_files / 2)
    ncols = 2 if n_bam_files > 1 else 1

    fig, axs = plt.subplots(nrows, ncols, figsize=(10 * ncols, 8 * nrows), squeeze=False)

    plotted_samples = []
    for idx, bam_file in enumerate(bam_files):
        sample = Path(bam_file).stem
        bam_file_path = os.path.join(bam_directory, bam_file)
        output_csv_path = os.path.join(output_csv, f"{sample}_noise.csv")

        average_depth = get_average_depth(bam_file_path, reference_name)
        if average_depth > 49:
            plotted_samples.append(bam_file)

    n_plotted_samples = len(plotted_samples)
    nrows = math.ceil(n_plotted_samples / 2)
    ncols = 2 if n_plotted_samples > 1 else 1

    fig, axs = plt.subplots(nrows, ncols, figsize=(10 * ncols, 8 * nrows), squeeze=False)

    min_base_proportion_threshold = 0.1

    plot_idx = 0
    for idx, bam_file in enumerate(bam_files):
        sample = Path(bam_file).stem
        bam_file_path = os.path.join(bam_directory, bam_file)
        output_csv_path = os.path.join(output_csv, f"{sample}_noise.csv")

        average_depth = get_average_depth(bam_file_path, reference_name)
        if average_depth > 49:
            row_idx = plot_idx // ncols
            col_idx = plot_idx % ncols
            ax = axs[row_idx, col_idx]

            save_noise_scores_to_csv(sample, bam_file_path, reference_name, output_csv_path, fig, ax, min_base_proportion=min_base_proportion_threshold)

            plot_idx += 1

    # Customize the remaining unused axes if there are any
    for idx in range(n_plotted_samples, nrows * ncols):
        row_idx = idx // ncols
        col_idx = idx % ncols
        ax = axs[row_idx, col_idx]
        ax.axis('off')

    customize_and_display_charts(fig, axs, plotted_samples)










'''
import pysam
import math
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
from matplotlib import rcParams

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
    df.to_csv(output_csv, index=False)

    # Set the Bodoni 72 font (assuming it's installed on your system)
    rcParams['font.family'] = 'Bodoni 72'
    rcParams['font.weight'] = 'bold'

    # Set up the Seaborn style to resemble Vox charts
    sns.set_style("whitegrid", {'grid.linestyle': '--'})
    sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})

    # Create a bar chart
    fig, ax = plt.subplots(figsize=(10, 8))
    bars = ax.bar(df['position'], df['noise'], color='#255F85', edgecolor='black', linewidth=1)

    # Add a small dot on top of every bar with a value higher than 0.5
    threshold = 0.5
    for idx, bar in enumerate(bars):
        bar_height = bar.get_height()
        if bar_height > threshold:
            ax.scatter(bar.get_x() + bar.get_width() / 2, bar_height, color='red', marker='o', s=30, zorder=2)


    # Customize chart appearance
    ax.set_title('Bar Chart', fontsize=18, color='#255F85')
    ax.set_xlabel('Position', fontsize=14, color='#255F85')
    ax.set_ylabel('Noise', fontsize=14, color='#255F85')
    ax.set_ylim(0, 1.0)


    # Set background color to light beige
    fig.patch.set_facecolor('#F5F5DC')
    ax.set_facecolor('#F5F5DC')

    # Remove default spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Change tick colors
    ax.tick_params(axis='x', colors='#255F85')
    ax.tick_params(axis='y', colors='#255F85')

    # Remove horizontal grid lines
    ax.yaxis.grid(False)

    # Add a thicker border around the chart area
    border = Rectangle((0, 0), 1, 1, transform=ax.transAxes,
                    fill=False, edgecolor='#255F85', linewidth=2.5)
    ax.add_patch(border)

    # Show the chart
    plt.show()
    #plt.savefig('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/plot.png', dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())









if __name__ == "__main__":
    sample = "25Neg_2"
    bam_file = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/bam/252301516_N.bam"
    reference_name = "A_HA_H1"
    output_csv = "/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/noise.csv"
    save_noise_scores_to_csv(sample, bam_file, reference_name, output_csv)

'''