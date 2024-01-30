# This script creates an animation of the strain against time plot, 
# eventually to be used with VisIt animations. 
#
# Note: Specific dimentions can be specified below, according to 
# the dimentions of your VisIt animation.

import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
import sys

# Check if the necessary arguments have been provided
if len(sys.argv) < 4:
    print("Usage: python script.py <output_directory> <save_video> <show plot>")
    sys.exit(1)

# Get the output directory and save video option from command line arguments
output_directory = sys.argv[1]
save_video = sys.argv[2].lower() == 'true'
show_plt = sys.argv[3].lower() == 'true'

# Get a list of all files that begin with "Rpsi4"
input_files = glob.glob("/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/Rpsi4*")

def valid_line(line):
    return not line.startswith("#")

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

    return data

# Load data from all files
data = [load_data(input_file) for input_file in input_files]

# Gather strain data
h_time = data[0][:, 0]
h_strain = data[3][:, 1]
strain_magnitude = np.real(h_strain)

# plot the entire line initially
fig, ax = plt.subplots(figsize=(772/80, 80/80))
line, = ax.plot(h_time, strain_magnitude, color='lightgray')
ani_line, = ax.plot([], [], color='blue')

ax.axis('off')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.set_xlim(h_time.min(), h_time.max())
ax.set_ylim(strain_magnitude.min(), strain_magnitude.max())

def update(i):
    ani_line.set_data(h_time[:i], strain_magnitude[:i])
    return ani_line,

# Create animation
ani = FuncAnimation(fig, update, frames=range(1, len(h_time)), blit=True)

# Save the animation
if save_video:
    Writer = writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    ani.save(os.path.join(output_directory, 'animation.mp4'), writer=writer)

print("Movie successfully saved to " + output_directory)

if show_plt:
    plt.show()