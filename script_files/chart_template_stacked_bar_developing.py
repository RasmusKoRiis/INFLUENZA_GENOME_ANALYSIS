import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib import rcParams

# Set the Bodoni 72 font (assuming it's installed on your system)
rcParams['font.family'] = 'Bodoni 72'
rcParams['font.weight'] = 'bold'

# Data
data = {
    'Sample': list(range(1, 11)),
    'HA': [25, 30, 20, 15, 35, 40, 30, 25, 10, 20],
    'NA': [15, 20, 10, 5, 25, 30, 20, 15, 5, 10],
    'MP': [10, 15, 5, 2, 20, 25, 15, 10, 2, 5],
    'NS': [5, 10, 2, 1, 15, 20, 10, 5, 1, 2],
    'NP': [20, 25, 15, 10, 30, 35, 25, 20, 5, 15],
    'PA': [10, 15, 5, 2, 20, 25, 15, 10, 2, 5],
    'PB1': [15, 20, 10, 5, 25, 30, 20, 15, 5, 10],
    'PB2': [10, 15, 5, 2, 20, 25, 15, 10, 2, 5]
}
df = pd.DataFrame(data)

# Set up the Seaborn style to resemble Vox charts
sns.set_style("whitegrid", {'grid.linestyle': '--'})
sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})

# Create a custom color palette with the provided colors
custom_palette = ['#FFC857', '#04A777', '#A41623', '#2D3047', '#255F85']

# Create a stacked bar chart
fig, ax = plt.subplots(figsize=(10, 8))
bars = []
for idx, fragment in enumerate(df.columns[1:]):
    if idx == 0:
        bars.append(ax.bar(df['Sample'], df[fragment], label=fragment, color=custom_palette[idx % len(custom_palette)]))
    else:
        bars.append(ax.bar(df['Sample'], df[fragment], bottom=df[df.columns[1:idx+1]].sum(axis=1),
                           label=fragment, color=custom_palette[idx % len(custom_palette)]))

# Customize chart appearance
ax.set_title('Stacked Bar Chart', fontsize=18, color='#255F85')
ax.set_xlabel('Samples', fontsize=14, color='#255F85')
ax.set_ylabel('Values', fontsize=14, color='#255F85')

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

# Move the legend outside the plot area, remove the background and border
legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=12)
frame = legend.get_frame()
frame.set_facecolor('none')
frame.set_edgecolor('none')

# Show the chart
plt.show()

