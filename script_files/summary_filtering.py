import pandas as pd
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep=",")

# Read the csv file into a DataFrame
df = pd.read_csv(input_file)

# Filter out the rows where the value in the Mapped_Reads column is less than 15
df = df[df['Mapped_Reads'] >= 15]

# Save the filtered DataFrame to a new csv file
df.to_csv(output_file, index=False)


# Save the final dataframe to a CSV file
df.to_csv(output_file, index=False)