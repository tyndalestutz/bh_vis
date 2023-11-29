import os, math, time, vtk
import numpy as np
from scipy.special import sph_harm, comb, lpmv

def old_sph_harm_array(l, m, s, radius_array, theta_array): #uses scipy.special sph_harm
    sph_harm_points = np.zeros((len(radius_array), len(theta_array)), dtype=np.complex128)
    phi = math.pi / 2  # phi is fixed to pi/2, why?
    spin_weight = (-1) ** s * math.sqrt((2 * l + 1) * math.factorial(l - m) / math.factorial(l + m) / (4 * math.pi)) # why?
    for j, radius in enumerate(radius_array):
        for i, theta in enumerate(theta_array):
            sph_harm_points[j, i] = spin_weight * scipy.special.sph_harm(m, l, theta, phi)
    return sph_harm_points

def spin_weight_sph_harm(l, m, s, theta, phi):
    sign = (-1)**(l + m - s)
    root_factor = math.sqrt((math.factorial(l + m) * math.factorial(l - m) * (2*l + 1)) / (4*math.pi * math.factorial(l + s) * math.factorial(l - s)))
    sine_and_phase = math.sin(theta/2) * math.exp(1j * m * phi)
    sum = 0
    for r in range(0, l - s + 1):
        sum += ((-1)**r) * comb(l - s, r) * comb(l + s, r + s - m) * ((math.cos(theta/2) / math.sin(theta/2))**(2*r + s - m))
    s_Y_lm = sign * root_factor * sine_and_phase * sum
    return s_Y_lm
def spin_weight_sph_harm2(l, m, s, theta, phi):
    sign = (-1)**(m)
    root_factor = math.sqrt(math.factorial(l + m) * math.factorial(l - m) * (2*l + 1)) / (4*math.pi * math.factorial(l + s) * math.factorial(l - s))
    phase = math.exp(1j * m * phi)
    legendre_term = lpmv(m, l, math.cos(theta))
    s_Y_lm = sign * root_factor * phase * legendre_term
    return s_Y_lm

def set_sph_harm_array(l, m, s, theta_array): #uses scipy.special sph_harm

    sph_harm_points = np.zeros((len(theta_array)), dtype=np.complex128)
    phi = math.pi / 2
    for i, theta in enumerate(theta_array):
        sph_harm_points[i] = spin_weight_sph_harm(l, m, s, phi, theta)
    return sph_harm_points

def main():
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
    theta_values = np.linspace(0, 2 * math.pi, numTheta, endpoint=False)
    our_old = old_sph_harm_array(l,m,s,radius_values,theta_values)
    new_1 = 1
    print(our_old)


if __name__ == "__main__":
    main()
    