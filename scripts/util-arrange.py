# The purpose of this script is to clean CSV data according to the time column. 
# This is a potentially vital step in preparing puncture tracker data.

import sys
import pandas as pd

# The first command-line argument is the input filename
input_file = sys.argv[1]

# Load the data from the CSV file
data = pd.read_csv(input_file)

# Sort the DataFrame by the 'time' column
sorted_data = data.sort_values('time')

# Generate the output filename based on the input filename
output_file = input_file.rsplit('.', 1)[0] + '_sorted.csv'

# Save the sorted DataFrame to the new CSV file
sorted_data.to_csv(output_file, index=False)