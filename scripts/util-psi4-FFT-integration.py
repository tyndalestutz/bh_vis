# Author: Zach Etienne
# Refactored by Nathan Nguyen for readability

import numpy as np
import os
import sys
from scipy.optimize import curve_fit


def read_psi4(file_name):
    """
    Reads an ASCII file with a header describing the real and imaginary parts of the data for each mode.
    Returns the data in a convenient format to access the real and imaginary parts given l, m values.

    Args:
        file_name (str): The name of the file to read.

    Returns:
        tuple: A tuple containing the time numpy array and a dictionary with keys (l, m) containing the data.
    """
    mode_data = {}
    time_data = []

    # Read file and ignore lines starting with #
    with open(file_name, "r") as file:
        lines = [line for line in file.readlines() if not line.startswith("#")]

    # Convert lines to arrays and sort by time
    data = np.array([list(map(np.float64, line.split())) for line in lines])
    data = data[np.argsort(data[:, 0])]

    # Remove duplicate times
    _, index = np.unique(data[:, 0], return_index=True)
    data = data[index]

    time_data = data[:, 0]

    # Extract l from the filename
    l = int(file_name.split("_")[1][1:].split("-")[0])

    # Loop through columns and store real and imaginary parts in mode_data
    for m in range(-l, l + 1):
        idx = 1 + 2 * (m + l)  # Calculate the index of the real part
        mode_data[(l, m)] = data[:, idx], data[:, idx + 1]

    return time_data, mode_data


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

    # Filter the data for t=100 to t=300
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
    if len(sys.argv) != 4:
        print("Usage: util-psi4-FFT-integration.py <path to gravitational wave Rpsi4 file> <path to output directory/> <l value>")
        sys.exit(1)

    file_name = sys.argv[1]
    output_directory = sys.argv[2]
    l_value = int(sys.argv[3])

    if not os.path.isfile(file_name):
        print(f"File {file_name} does not exist. Ending...")
        return

    # Read the data from the file
    time, mode_data = read_psi4(file_name)

    # Loop over each m value
    for m in range(-l_value, l_value + 1):
        real, imag = mode_data[(l_value, m)]
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

        output_file = output_directory + "Rpsi4_l_" + str(l_value) + "_m_" + str(m) + "-r0100.0.txt"
        with open(output_file, "w") as file:
            file.write("# Time    Second_Integral_Real    Second_Integral_Imag\n")
            for t, real, imag in zip(time, intg_dt2_real, intg_dt2_imag):
                file.write(f"{t:.15f} {real:.15f} {imag:.15f}\n")

        print(f"Second time integral data has been saved to {output_file}")


if __name__ == "__main__":
    main()
    