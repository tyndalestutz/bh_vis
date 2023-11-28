import os, math, time
import re
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from scipy.special import sph_harm, factorial, comb, lpmv
from scipy.optimize import curve_fit
import spherical
import quaternionic

input_directory = r'/home/guest/Documents/Users/Seth/bh_vis/scripts/data'


# Plot parameters
num_radius_pts = 450
num_theta_pts = 180
display_radius = 300
R_ext = 100
ell_max = 8
l_vals = range(2,ell_max+1)
m = 2
s = -2
scale_factor = 600
omit_radius = 4
plot_strain = True
perform_integration = True


# Find Wigner matrix for every l from 1 to ell_max
wigner = spherical.Wigner(ell_max)



# This function is currently the most accurate with the Mathematica notebook I tested:
# https://www.black-holes.org/SpinWeightedSphericalHarmonics.nb
def spin_weight_sph_harm(l, m, s, theta, phi):
    """
    Finds spin weighted spherical harmonics for a given l, m, s, and array of theta values.
    Based off the Goldberg et al. formula found on Wikipedia.

    Theta is polar angle (0 to pi), phi is azimuthal (0 to 2pi)

    Returns a single complex value based on parameters.
    """
    sign = (-1)**(l + m - s)
    root_factor = (factorial(l + m) * factorial(l - m) * (2*l + 1)) / (4*np.pi * factorial(l + s) * factorial(l - s))
    root_factor = np.sqrt(root_factor)
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
    root_factor = (factorial(l + m) * factorial(l - m) * (2*l + 1)) / (4*np.pi * factorial(l + s) * factorial(l - s))
    root_factor = np.sqrt(root_factor)
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

r_vals = np.linspace(0, display_radius, num_radius_pts)
theta_vals = np.linspace(0, 2*np.pi, num_theta_pts, endpoint=False)
#sph_array = set_sph_harm_array(2, 2, -2, theta_vals)

# sort psi4 modes and store into respective arrays
psi4_file_pattern = r"Rpsi4_l([\d.]+)-r0100.0.txt"
# organized data will be stored in this dictionary
psi4_data = {}
for current_file in os.listdir(input_directory):
    match = re.search(psi4_file_pattern, current_file)
    if match:
        l = int(match.group(1))
        file_path = os.path.join(input_directory, current_file)
        with open(file_path, 'r') as file:
            lines = [line for line in file.readlines() if not line.startswith('#')]
            current_data = np.array([list(map(float, line.split())) for line in lines]) # 2d, [row, column]
            data_times = current_data[:, 0] # store first column for psi4 times
            data_times_reshaped = data_times.reshape(-1, 1)
            for m in range(-l, l + 1):
                psi4_real = current_data[:, 2*(m + l) + 1]
                psi4_imag = current_data[:, 2*(m + l) + 2]
                print(psi4_imag)
                #####psi4_complex = np.vectorize(complex)(psi4_real, psi4_imag)
                # the current l, m mode psi4 data is now stored and added to dictionary
                #####combined_data = np.column_stack((data_times, psi4_complex))
                combined_data = np.column_stack((data_times, psi4_real, psi4_imag))
                # each dictionary item is a 2d array with each row: [time, psi4 (complex)]
                psi4_data[f'l{l}_m{m}'] = combined_data

# loop through all organized psi4 modes and multiply by corresponding
# spin-weighted spherical harmonic
# strain_value = S.real * h_tR.real - S.imag * h_tR.imag
for theta_pt, theta_val in enumerate(theta_vals):
    #####psi4_data[f's-imposed_theta{theta_pt}'] = np.zeros((len(psi4_data[f'l{l_vals[0]}_m0'][:, 0]), 2), dtype=complex)
    psi4_data[f's-imposed_theta{theta_pt}'] = np.zeros((len(psi4_data[f'l{l_vals[0]}_m0'][:, 0]), 2), dtype=np.complex128)
    psi4_data[f's-imposed_theta{theta_pt}'][:, 0] = psi4_data[f'l{l_vals[0]}_m0'][:, 0]
    psi4_data[f'strain_full_theta{theta_pt}'] = np.zeros((len(psi4_data[f'l{l_vals[0]}_m0'][:, 0]), 3))

'''
NOTE: uses old version of swsh
for l in l_vals:
    for m in range(-l, l + 1):
        psi4_time = psi4_data[f'l{l}_m{m}'][:, 0]
        sph_array = set_sph_harm_array(l, m, s, theta_vals)
        for theta_pt, theta_val in enumerate(theta_vals):
            s_factor = sph_array[theta_pt]
            psi4_real = psi4_data[f'l{l}_m{m}'][:, 1]
            psi4_imag = psi4_data[f'l{l}_m{m}'][:, 2]
            psi4_with_sph_harm = (psi4_real * s_factor.real - psi4_imag * s_factor.imag) + 1j * ((psi4_imag * s_factor.real) + (psi4_real * s_factor.imag))
            ###psi4_with_sph_harm = [x.real * s_factor.real -
                                  ###x.imag * s_factor.imag
                                  ###for x in psi4_data[f'l{l}_m{m}'][:, 1]]
            psi4_data[f'l{l}_m{m}_theta{theta_pt}'] = np.column_stack((psi4_time, psi4_with_sph_harm))
            psi4_data[f's-imposed_theta{theta_pt}'][:, 1] += psi4_data[f'l{l}_m{m}_theta{theta_pt}'][:, 1]
            if theta_pt == 0:
                print(f'added mode l={l}, m={m}')
'''

# polar angle
phi = np.pi/2
#phi = 90
#theta_vals_deg = np.degrees(theta_vals)
for theta_pt, theta_val in enumerate(theta_vals):
    R = quaternionic.array.from_spherical_coordinates(theta_val, phi)
    #spin weighted spherical harmonics array
    Y = wigner.sYlm(s, R)
    for l in l_vals:
        for m in range(-l, l + 1):
            psi4_time = psi4_data[f'l{l}_m{m}'][:, 0]
            swsh_index = l**2 + m + l # this index pulls from swsh array based off of current l, m
            #s_factor = Y[swsh_index]
            s_factor = Y[wigner.Yindex(l, m)]
            psi4_real = psi4_data[f'l{l}_m{m}'][:, 1]
            psi4_imag = psi4_data[f'l{l}_m{m}'][:, 2]
            psi4_with_swsh = (psi4_real * s_factor.real - psi4_imag * s_factor.imag) + 1j * ((psi4_imag * s_factor.real) + (psi4_real * s_factor.imag))
            
            psi4_data[f'l{l}_m{m}_theta{theta_pt}'] = np.column_stack((psi4_time, psi4_with_swsh))
            psi4_data[f's-imposed_theta{theta_pt}'][:, 1] += psi4_data[f'l{l}_m{m}_theta{theta_pt}'][:, 1]
            if theta_pt == 0:
                print(f'added mode l={l}, m={m}')


# Some notes to keep track of how swsh affects amplitudes
'''
l2m2: beginning peaks at 0.00015, ending peaks at 0.0035 - spherical harmonics reduces by a factor of 20
l2m0: beginning peaks at 0.0001, ending peaks at 0.0002 - spherical harmonics reduces by a factor of 5
l4m0: beginning peaks at 0.0006, ending peaks at - spherical harmonics reduces by a factor of 3.6
l6m0: beginning peaks at 0.0004, ending peaks at - spherical harmonics reduces by a factor of 3.1
l8m0: beginning peaks at 0.0002, ending peaks at - spherical harmonics reduces by a factor of 2.7
adding up all the modes should have a beginning peak no greater than 0.04 or so
'''

print('Psi4 organizaion and spin-weighted spherical harmonic\n muliplication finished. Integrating...')


'''
Next: integrate n times for our n superimposed psi4 wave modes, corresponding to our n
theta values varying in spin weighted spherical harmonics.

'''
def compute_first_derivative(time, data):
    """
    Calculates the time derivative of the input data using a second-order finite difference stencil.

    Args:
        time (numpy.ndarray): A numpy array containing time values.
        data (numpy.ndarray): A numpy array containing the data to be differentiated.

    Returns:
        numpy.ndarray: A numpy array containing the time derivative of the input data.
    """
    dt = time[1] - time[0]
    derivative = np.zeros_like(data)
    # Second-order in the interior
    derivative[1:-1] = (data[2:] - data[:-2]) / (2 * dt)
    # Drop to first-order at the endpoints
    derivative[0] = (data[1] - data[0]) / dt
    derivative[-1] = (data[-1] - data[-2]) / dt

    return derivative

def process_wave_data(time, real, imag):
    """
    Calculates the cumulative phase and amplitude of a gravitational wave signal.

    Args:
        time (numpy.ndarray): A numpy array containing time values.
        real (numpy.ndarray): A numpy array containing the real part of the signal.
        imag (numpy.ndarray): A numpy array containing the imaginary part of the signal.

    Returns:
        tuple: A tuple containing three numpy arrays (amplitude, cumulative phase, phase time-derivative).
    """
    amplitudes = np.sqrt((real**2 + imag**2), dtype=np.float64)
    phases = np.arctan2(imag, real, dtype=np.float64)

    cycles = 0
    cumulative_phases = np.zeros_like(phases)
    last_phase = phases[0]

    for index, phase in enumerate(phases):
        #Identify phase wrapping
        if phase - last_phase > np.pi:
            cycles -= 1
        if phase - last_phase < -1 * np.pi:
            cycles += 1

        cumulative_phases[index] = phase + cycles * 2 * np.pi
        last_phase = phase

    cumulative_phase_derivatives = compute_first_derivative(time, cumulative_phases)

    return  amplitudes, cumulative_phases, cumulative_phase_derivatives

def fit_quadratic_and_output_min_omega(time, omega):
    def quadratic(x, a, b, c):
        return a * x**2 + b * x + c

    # Filter the data for t=200 to t=400. Note: non-pysical data at t=0 to t=120
    time_filtered = time[(time >= 100) & (time <= 300)]
    omega_filtered = omega[(time >= 100) & (time <= 300)]

    # Fit a quadratic curve to the omega data using nonlinear least squares
    params, _ = curve_fit(quadratic, time_filtered, omega_filtered)
    a, b, c = params

    # Find the extremum value of the quadratic curve and estimate omega(0)
    extremum_x = -b / (2 * a)
    omega_min_quad_fit = np.abs(quadratic(extremum_x, a, b, c))
    omega_at_time_zero = np.abs(quadratic(0.0, a, b, c))

    print(f"The extremum of the quadratic curve occurs at t = {extremum_x:.15f} with omega = {omega_min_quad_fit:.15f} . implied omega(t=0) = {omega_at_time_zero:.15f}")
    return np.abs(omega_at_time_zero)

def perform_complex_fft(time, real, imag):
    """
    Performs a complex Fast Fourier Transform (FFT) on the input time, real, and imaginary data.

    Args:
        time (numpy.ndarray): A numpy array containing time values.
        real (numpy.ndarray): A numpy array containing the real part of the signal.
        imag (numpy.ndarray): A numpy array containing the imaginary part of the signal.

    Returns:
        tuple: A tuple containing two numpy arrays (angular_frequencies, fft_data).
    """
    # Combine the real and imaginary data into a single complex signal
    complex_signal = real + 1j * imag
    
    # Perform the complex FFT
    fft_data = np.fft.fft(complex_signal)

    # Calculate the frequency values
    dt = time[1] - time[0]
    n = len(time)
    angular_frequencies = np.fft.fftfreq(n, d=dt) * 2 * np.pi

    return angular_frequencies, fft_data


def main():
    """
    Main function that reads the gravitational wave data file, processes the data,
    and saves the output to a file. The input filename is provided via command line.
    """

    output_directory = '/home/guest/Documents/Users/Seth/bh_vis/scripts/r0100.0/S-imposed_strain'

    # NOTE: instead of reading file we pull from psi4 dictionary

    # Loop over each theta value to integrate
    for theta_pt, theta_val in enumerate(theta_vals):
        time = psi4_data[f's-imposed_theta{theta_pt}'][:, 0].real
        real = psi4_data[f's-imposed_theta{theta_pt}'][:, 1].real
        imag = psi4_data[f's-imposed_theta{theta_pt}'][:, 1].imag
        _, _, phase_dt = process_wave_data(time, real, imag)
        omega_0 = fit_quadratic_and_output_min_omega(time, phase_dt)
        # Perform FFT
        fft_ang_freqs, fft_result = perform_complex_fft(time, real, imag)

        # Fixed Frequency Integration
        # Just below Eq. 27 in https://arxiv.org/abs/1006.1632
        for i, omega in enumerate(fft_ang_freqs):
            if np.abs(omega) <= omega_0:
                fft_result[i] *= 1 / (1j * omega_0) ** 2
            else:
                fft_result[i] *= 1 / (1j * np.abs(omega)) ** 2

        # Now perform the inverse FFT
        intg_dt2_complex = np.fft.ifft(fft_result)

        # Separate the real and imaginary parts of the second time integral
        intg_dt2_real = np.real(intg_dt2_complex)
        intg_dt2_imag = np.imag(intg_dt2_complex)



        output_file = output_directory + f"/S-imposed_strain_theta{theta_pt}-r0100.0.txt"
        with open(output_file, "w") as file:
            file.write("# Time    Second_Integral_Real    Second_Integral_Imag\n")
            for t, real, imag in zip(time, intg_dt2_real, intg_dt2_imag):
                file.write(f"{t:.15f} {real:.15f} {imag:.15f}\n")

        psi4_data[f'strain_full_theta{theta_pt}'][:, 0] = psi4_time
        psi4_data[f'strain_full_theta{theta_pt}'][:, 1] = intg_dt2_real
        psi4_data[f'strain_full_theta{theta_pt}'][:, 2] = intg_dt2_imag



        print(f"Second time integral data has been saved to {output_file}")


if __name__ == "__main__" and perform_integration:
    main()

# Temporary 2d animations in matplotlib for testing
'''
fig, ax = plt.subplots()
def update(frame):
    line.set_ydata(psi4_data[f'strain_full_theta{frame}'][:, 1].real)
    ax.set_title(f'theta point {frame}')
    return line,
line, = ax.plot(psi4_data['strain_full_theta0'][:, 0], psi4_data['strain_full_theta0'][:, 1].real)
animation = FuncAnimation(fig, update, frames=range(180), interval=5)


#plt.plot(psi4_data['s-imposed_theta0'][:, 0].real, psi4_data['s-imposed_theta0'][:, 1].real, color = 'b', alpha=0.5)
#plt.plot(psi4_data['l2_m2_theta0'][:, 0].real, psi4_data['l2_m2_theta0'][:, 1].real, color='r', alpha=0.5)
#plt.plot(psi4_data['l2_m2'][:, 0], psi4_data['l2_m2'][:, 1].real, color='b', alpha=0.5)
plt.grid(True)
plt.show()
'''
