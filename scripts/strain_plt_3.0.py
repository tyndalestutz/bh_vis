# Author: Tyndale Stutzman
# Last Modified: 07.06.23
#
# Note: This script is currently functional. It's purpose is to generate a bare bones, basic plot
#   for each mode and save it as a png. 

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

for data, l, m in data_modes:
    h_time = data[:, 0]
    h_strain = data[:, 1]
    strain_magnitude = np.real(h_strain)

    # Create the plot
    fig, ax = plt.subplots()

    # Plot the data
    ax.plot(h_time, strain_magnitude)

    # Set titles for the axes
    ax.set_xlabel('Time')
    ax.set_ylabel('Strain')

    # Set title for the plot
    ax.set_title(f'Mode (l, m): ({l}, {m})')

    # Save the figure
    #plt.savefig(f"{output_directory}strain_mode_l{l}_m{m}.png")
    
    #plt.close(fig)
