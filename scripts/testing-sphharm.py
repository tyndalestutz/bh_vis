import os, math, time, vtk
import numpy as np
from scipy.special import sph_harm, factorial, comb, lpmv


# Plot parameters
numRadius = 450  # Number of points in the radial direction
numTheta = 180  # Number of points per circle
display_radius = 300  # Mesh radius
R_ext = 100  # Extraction radius
l = 2
m = 2
s = -2
scale_factor = 600
omit_radius = 4


radius_values = np.linspace(0, display_radius, numRadius)
theta_values = np.linspace(0, 2 * np.pi, numTheta, endpoint=False)


def old_sph_harm_array(l, m, s, radius_array, theta_array): #uses scipy.special sph_harm
    sph_harm_points = np.zeros((len(radius_array), len(theta_array)), dtype=np.complex128)
    phi = np.pi / 2  # phi is fixed to pi/2, why?
    spin_weight = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * math.factorial(l - m) / math.factorial(l + m)) # why?
    for j, radius in enumerate(radius_array):
        for i, theta in enumerate(theta_array):
            sph_harm_points[j, i] = spin_weight * sph_harm(m, l, theta, phi)
    return sph_harm_points

def spin_weight_sph_harm(l, m, s, theta, phi):
    """
    Finds spin weighted spherical harmonics for a given l, m, s, and array of theta values.
    Based off the Goldberg et al. formula found on Wikipedia.

    Theta is polar angle (0 to pi), phi is azimuthal (0 to 2pi)

    Returns a single complex value based on parameters.
    """
    sign = (-1)**(l + m - s)
    root_factor = np.sqrt((factorial(l + m) * factorial(l - m) * (2*l + 1)) / (4*np.pi * factorial(l + s) * factorial(l - s)))
    sine_and_phase = np.sin(theta/2) * np.exp(1j * m * phi)
    sum = 0
    for r in range(0, l - s + 1):
        sum += ((-1)**r) * comb(l - s, r) * comb(l + s, r + s - m) * ((np.cos(theta/2) / np.sin(theta/2))**(2*r + s - m))
    s_Y_lm = sign * root_factor * sine_and_phase * sum
    return s_Y_lm

def spin_weight_sph_harm2(l, m, s, theta, phi):
    """
    Finds spin weighted spherical harmonics for a given l, m, s, and array of theta values.
    Based off the Goldberg et al. formula found on Wikipedia.

    Theta is polar angle (0 to pi), phi is azimuthal (0 to 2pi)

    Returns a single complex value based on parameters.
    """
    sign = (-1)**(m)
    root_factor = np.sqrt(factorial(l + m) * factorial(l - m) * (2*l + 1)) / (4*np.pi * factorial(l + s) * factorial(l - s))
    phase = np.exp(1j * m * phi)
    legendre_term = lpmv(m, l, np.cos(theta))
    s_Y_lm = sign * root_factor * phase * legendre_term
    return s_Y_lm

# here theta is azimuthal angle (0 to 2pi) and phi is the polar angle (0 to pi) (keep at pi/2 for 2d representation)
def set_sph_harm_array(l, m, s, theta_array): #uses scipy.special sph_harm
    """
    Finds spin-weighted spherical harmonics for a given l, m, s, and array of theta values. 
    Uses scipy.special sph_harm for initial spherical harmonic factor.

    Theta is our azimuthal angle going from 0 to 2pi,
    for now our polar angle phi is fixed to pi/2 for 2d representation.

    Returns a 1d array of spherical harmonics for varying theta.
    """
    sph_harm_points = np.zeros((len(theta_array)), dtype=np.complex128)
    phi = np.pi / 2
    #spin_weight = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * math.factorial(l - m) / math.factorial(l + m))
    for i, theta in enumerate(theta_array):
        #spin_weight = np.exp(1j * m * theta)
        #spin_weight = ((-1)**(l + m - s))*np.sqrt((factorial(l + m) * factorial(l + m)) / (factorial(l + s) * factorial(l - s)))
        # NOTE: This function switches phi and theta definitions from spin_weight_sph_harm()
        sph_harm_points[i] = spin_weight_sph_harm(l, m, s, phi, theta)
        #sph_harm_points[i] = spin_weight * sph_harm(m, l, theta, phi) 

    return sph_harm_points

our_old = old_sph_harm_array(l,m,s,radius_values,theta_values)
new_1 = 