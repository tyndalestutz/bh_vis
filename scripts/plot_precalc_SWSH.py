import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle

# Load the dictionary
with open('/home/guest/Documents/Users/Tyndale/bh_repos/bhah_waveform_analysis/r0100.0/many_modes/test_calc.pkl', 'rb') as f:
    sph_harm_dict = pickle.load(f)

# Now you can use it as a dictionary
print(sph_harm_dict.keys())
# print(sph_harm_dict[(2, 2)])

# Parameters
numTheta = 180
numRadius = 200
display_radius = 100

def plot_sph_harm(l, m, sph_harm_dict, display_radius):
    # Fetch the spherical harmonic array from the dictionary
    sph_harm_points = sph_harm_dict[(l, m)]
    
    # Create the meshgrid
    theta_vals = np.linspace(0, 2*np.pi, sph_harm_points.shape[0])
    radius_vals = np.linspace(0, display_radius, sph_harm_points.shape[1])
    theta, radius = np.meshgrid(theta_vals, radius_vals)
    x = radius * np.sin(theta) 
    y = radius * np.cos(theta)
    
    # Prepare the real part of the array, ensuring the shape matches the other arrays
    z_real = np.real(sph_harm_points.T)  # Transpose the array

    # Print the shapes to check them
    print(f"x shape: {x.shape}")
    print(f"y shape: {y.shape}")
    print(f"z_real shape: {z_real.shape}")

    # Plot the real part
    fig = plt.figure(figsize=(12, 6))

    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, z_real, rstride=1, cstride=1, color='c', alpha=0.6, linewidth=0)
    ax.set_title(f'Real part of Spherical Harmonic: l={l}, m={m}')

    plt.show()

# To plot spherical harmonics of a specific mode:
l, m = 2, 2  # example values
plot_sph_harm(l, m, sph_harm_dict, display_radius)
