import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers

# File parameters
folder_name = "sum_test"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/"

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

# Check that all time series match
#for i in range(len(data) - 1):
    #assert np.array_equal(data[i][0], data[i+1][0]), "Time series do not match"

# Sum up strain data
h_time = data[0][:, 0]
#h_strain = sum(item[:, 2] for item in data)
h_strain = data[3][:, 1]
strain_magnitude = np.real(h_strain)

###

# Create the plot
fig, ax = plt.subplots()

# Plot the data
ax.plot(h_time, strain_magnitude)

# Set titles for the axes
ax.set_xlabel('Time')
ax.set_ylabel('Strain')

###

'''
# plot the entire line initially
line, = ax.plot(h_time, strain_magnitude, color='lightgray')
ani_line, = ax.plot([], [], color='blue')
fig, ax = plt.subplots(figsize=(772/80, 80/80))

# Turn off axes
ax.axis('off')

# Remove bounding box
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.spines['left'].set_visible(False)

# Set the axes limits to maximize the plot size
ax.set_xlim(h_time.min(), h_time.max())
ax.set_ylim(strain_magnitude.min(), strain_magnitude.max())

def update(i):
    ani_line.set_data(h_time[:i], strain_magnitude[:i])
    return ani_line,

# Create animation
ani = FuncAnimation(fig, update, frames=range(1, len(h_time)), blit=True)

# Save the animation
#Writer = writers['ffmpeg']
#writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
#ani.save('animation.mp4', writer=writer)
'''
plt.show()