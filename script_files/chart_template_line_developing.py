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
    'Categories': [1, 2, 3, 4],
    'Series 1': [25, 45, 30, 15],
    'Series 2': [15, 35, 25, 10],
}
df = pd.DataFrame(data)

# Set up the Seaborn style to resemble Vox charts
sns.set_style("whitegrid", {'grid.linestyle': '--'})
sns.set_context("notebook", font_scale=1.25, rc={"lines.linewidth": 2.5})

# Create a custom color palette with the provided colors
custom_palette = ['#FFC857', '#04A777', '#A41623', '#2D3047', '#255F85']

# Create a line chart
fig, ax = plt.subplots(figsize=(8, 8))  # Set the plot size to be square
for idx, series in enumerate(df.columns[1:]):
    sns.lineplot(data=df, x='Categories', y=series, ax=ax, linewidth=2.5, color=custom_palette[idx], label=series)

# Customize chart appearance
ax.set_title('Vox-Style Line Chart', fontsize=18, color='#255F85')
ax.set_xlabel('Categories', fontsize=14, color='#255F85')
ax.set_ylabel('Values', fontsize=14, color='#255F85')

# Set background color to light beige
fig.patch.set_facecolor('#F5F5DC')  # Light beige color (you can modify this to your preferred color)
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
