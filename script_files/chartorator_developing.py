
import pandas as pd
import plotly.express as px

# Read the data into a pandas DataFrame
df = pd.read_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/mutation/A_H1_MP_primer_alignment.csv')


# Create the heatmap using Plotly Express
fig = px.imshow(df.pivot_table(values='Total_Mis', index='Primer', columns='SequenceID', aggfunc='first'),
                labels=dict(x='SequenceID', y='Primer', color='Total_Mis'),
                x=df['SequenceID'].unique(), y=df['Primer'].unique(),text_auto=True, )

# Set the title and axis labels
fig.update_layout(
    title='Heatmap of Total_Mis',
    xaxis_title='SequenceID',
    yaxis_title='Primer',
    xaxis=dict(showgrid=True, gridcolor='rgba(200, 200, 200, 1)'),  # Add x-axis grid
    yaxis=dict(showgrid=True, gridcolor='rgba(200, 200, 200, 1)')   # Add y-axis grid
)

fig.show()

# Save the figure as an HTML file
fig.write_html('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results_dev/heatmap.html')

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_csv('/Users/rasmuskopperudriis/Coding/work/new_influenza_pipline/results/new_influenza_pipline_final.csv')


# Assuming the DataFrame named df is already loaded with your CSV data

# Extract the mutation columns
mutation_columns = [col for col in df.columns if "Mutations" in col]

# Initialize an empty dictionary to store mutation counts
mutation_counts = {}

# Iterate over mutation columns and count occurrences of each mutation
for col in mutation_columns:
    for mutation_list in df[col].dropna():
        mutations = mutation_list.split(';')
        for mutation in mutations:
            if mutation not in mutation_counts:
                mutation_counts[mutation] = 1
            else:
                mutation_counts[mutation] += 1

# Convert mutation counts into a DataFrame
mutations_df = pd.DataFrame(list(mutation_counts.items()), columns=['Mutation', 'Count'])

# Create a pivot table with mutation columns and their counts
pivot_table = pd.pivot_table(mutations_df, values='Count', index='Mutation', columns=['Mutation'], fill_value=0)

# Create a heatmap using Seaborn
plt.figure(figsize=(12, 12))
sns.set(font_scale=0.8)
sns.heatmap(pivot_table, cmap='coolwarm', linewidths=0.5, annot=True, fmt="d", cbar=False)

# Customize the plot
plt.title("Mutation Frequencies")
plt.xlabel("Mutations")
plt.ylabel
plt.show()