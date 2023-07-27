# Created with GPT4, this script sums of all modes, excluding (l, 0) modes, 
# and compares them against the control mode, (2, 2).

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import re

# File parameters
folder_name = "sum_test"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/"

# Get a list of all files that begin with "Rpsi4"
input_files = glob.glob("/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/Rpsi4*")

def valid_line(line):
    return not line.startswith("#")

def get_mode(filename):
    match = re.search(r"_l_(\d+)_m_(-?\d+)", filename)
    if match:
        l_value, m_value = map(int, match.groups())
        return l_value, m_value
    else:
        return None, None

def load_data(file_name):
    with open(file_name, 'r') as f:
        valid_lines = [line for line in f if valid_line(line)]
        data = [list(map(float, line.split())) for line in valid_lines]

    # Determine expected number of columns based on the first row
    num_columns = len(data[0])

    # Filter out rows with a different number of columns
    data = [row for row in data if len(row) == num_columns]

    data = np.array(data)
    data = np.unique(data, axis=0)

    l_value, m_value = get_mode(file_name)
    return data, l_value, m_value

# Load data from all files
data_modes = [load_data(input_file) for input_file in input_files]

# Sort data_modes by l, then m
data_modes.sort(key=lambda x: (x[1], x[2]))

# Initialize sum_strain
sum_strain = np.zeros_like(data_modes[0][0][:, 1])

for data, l, m in data_modes:
    if m != 0:  # Exclude modes where m = 0
        sum_strain += np.real(data[:, 1])

h_time = data_modes[0][0][:, 0]

# Define specific mode to plot
l_specific = 2
m_specific = 2

# Find data for specific mode
specific_mode_data = None
for data, l, m in data_modes:
    if l == l_specific and m == m_specific:
        specific_mode_data = np.real(data[:, 1])
        break

# Create the plot
fig, ax = plt.subplots()

# Plot the sum of strain data
ax.plot(h_time, sum_strain, label='Sum of Strain excluding (l, 0) modes')

# Plot the strain data for the specific mode, if found
if specific_mode_data is not None:
    ax.plot(h_time, specific_mode_data, label=f'Control mode (l, m): ({l_specific}, {m_specific})', color='red')

# Set titles for the axes
ax.set_xlabel('Time')
ax.set_ylabel('Strain')
ax.legend()
plt.show()