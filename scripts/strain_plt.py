
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers

# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test10_polar_zeroR_r=300/"

# Reads inputted strain data - returns data file row length, initial strain data (complex), spin-weighted spherical harmonics (complex), and time values
def initialize():
    if os.path.exists(output_directory):
        if len(os.listdir(output_directory)) != 0:
            answer = input(f"Data already exists at {output_directory}. Overwrite it? You cannot undo this action. (Y/N) ")
            if answer.capitalize() == "Y":
                for file in os.listdir(output_directory):
                    os.remove(f"{output_directory}/{file}")
            else:
                print("Exiting Program. Change output directory to an empty directory.")
                exit()
    else:
        last_slash_index = output_directory.rfind("/")
        super_directory = output_directory[:last_slash_index] if last_slash_index != -1 else output_directory
        if os.path.exists(super_directory):
            os.makedirs(output_directory)
        else:
            print(f"Error: {super_directory} does not exist.")
            exit()

    def valid_line(line):
        return not line.startswith("#")

    with open(input_file, 'r') as f:
        strain_data = np.array([list(map(float, line.split())) for line in f if valid_line(line)])

    strain_data = np.unique(strain_data, axis=0)
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j * h_imag

    length = len(h_strain)

    return length, h_strain, h_time

length, h_strain, h_time = initialize()

######################

fig, ax = plt.subplots(figsize=(772/80, 80/80))

strain_magnitude = np.real(h_strain)

# plot the entire line initially
line, = ax.plot(h_time, strain_magnitude, color='lightgrey')
ani_line, = ax.plot([], [], color='blue')

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
Writer = writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
ani.save('animation.mp4', writer=writer)

plt.show()
##########################