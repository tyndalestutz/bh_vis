import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm
from mpl_toolkits.mplot3d import Axes3D

# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test10_polar_zeroR_r=300/"

# Parameters
numTheta = 180
numRadius = 200
display_radius = 100 
R_ext = 100 

def set_sph_harm_array(l_min, l_max, numTheta, numRadius, display_radius):
    # Initialize the spherical harmonic dictionary
    sph_harm_dict = {}

    # Loop over all l values
    for l in range(l_min, l_max + 1):
        # Loop over all m values
        for m in range(-l, l + 1):
            # Initialize the spherical harmonic array for this m
            sph_harm_points = np.zeros((numTheta, numRadius), dtype=np.complex128)

            # Loop over all points in the grid
            for j in range(numRadius):
                for i in range(numTheta):
                    # Calculate the radial and angular coordinates
                    radius = display_radius * j / (numRadius - 1)
                    theta = 2 * np.pi * i / numTheta
                    
                    # Calculate the spherical coordinates
                    r = radius
                    theta = theta
                    phi = np.pi / 2  # phi is fixed to pi/2, similar to the original script

                    # Calculate the spherical harmonic
                    Y_lm = sph_harm(m, l, theta, phi)

                    # Adjust for spin-weighted
                    s = -2  # Assuming s=-2 as it's common in gravitational wave analysis
                    Y = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * math.factorial(l - m) / math.factorial(l + m)) * Y_lm

                    # Store the spherical harmonic in the array
                    sph_harm_points[i, j] = Y

            # Store the array in the dictionary with key (l, m)
            sph_harm_dict[(l, m)] = sph_harm_points

    return sph_harm_dict

sph_harm_dict = set_sph_harm_array(2, 8, numTheta, numRadius, display_radius)

import pickle

# Save the dictionary
with open('/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/test_calc.pkl', 'wb') as f:
    pickle.dump(sph_harm_dict, f)