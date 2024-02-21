import os
from typing import List, Tuple
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

OUTPUT_DIR = "visit/output"
INPUT_DIR = "r100"
FILE_PATTERN = "_l[MODE=L]-"
ELL_MIN = 2
ELL_MAX = 8
EXT_RAD = 100
INTERVAL = 200
CUTTOFF_FACTOR = 0.75


def psi4_fft_to_strain() -> np.ndarray:
    """
    Calculates the strain modes from PSI4 data using FFT.

    Returns:
        A numpy array of numpy arrays representing the strain modes.

    Raises:
        IOError: If there is an error reading the PSI4 data.
        ValueError: If the lengths of the time and data arrays are not equal.
    """

    try:
        time_arr, psi4_modes_data = read_psi4_dir()
    except IOError as e:
        raise IOError(f"Error reading PSI4 data: {e}") from e

    # Get minimum frequency cutoff from l=m=2 mode
    min_omega_l2m2 = extract_min_omega_ell2_m2(
        time_arr, psi4_modes_data[get_mode_index(2, 2)]
    )
    freq_cutoff = CUTTOFF_FACTOR * min_omega_l2m2

    # Initialize arrays for strain modes and their second time derivatives
    strain_modes = np.zeros_like(psi4_modes_data, dtype=np.complex128)
    strain_modes_ddot = np.zeros_like(psi4_modes_data, dtype=np.complex128)

    # Calculate frequency list for FFT
    freq_list = np.fft.fftfreq(len(time_arr), time_arr[1] - time_arr[0]) * 2 * np.pi

    # Loop over modes and calculate strain
    mode_idx = 0
    for l in range(ELL_MIN, ELL_MAX + 1):
        for m in range(-l, l + 1):
            mode_data = psi4_modes_data[get_mode_index(l, m)]

            # Apply FFT and filter, see Eq. 27 in https://arxiv.org/abs/1006.1632
            fft_result = np.fft.fft(mode_data)
            for i, omega in enumerate(freq_list):
                if np.fabs(omega) <= np.fabs(freq_cutoff):
                    fft_result[i] *= 1 / (1j * freq_cutoff) ** 2
                else:
                    fft_result[i] *= 1 / (1j * omega) ** 2

            # Inverse FFT to get strain
            strain_modes[mode_idx] = np.fft.ifft(fft_result)

            # Calculate second time derivative
            strain_modes_ddot[mode_idx] = second_time_derivative(
                time_arr, strain_modes[mode_idx]
            )
            #plt.plot(time_arr, mode_data[mode_idx] - strain_modes_ddot[mode_idx])

            mode_idx += 1
        #plt.show()
    return strain_modes


def read_psi4_dir() -> tuple[np.ndarray, np.ndarray]:
    """
    Reads data from psi4 output directory and returns time and mode data.

    Returns:
        tuple[np.ndarray, np.ndarray]:
            - time_data: Array of time values (shape: (n_times,)).
            - mode_data: Array of complex modes (shape: (n_modes, n_times, 2*l+1)).
    """

    time_data: list[np.ndarray] = []
    psi4_modes_data: list[np.ndarray] = []

    for l in range(ELL_MIN, ELL_MAX + 1):
        filepath = find_file_for_l(l)
        with open(filepath, "r", encoding="utf-8") as file:
            lines = [line for line in file.readlines() if not line.startswith("#")]
        data = np.array([list(map(np.float64, line.split())) for line in lines])

        # Extract time data once for all modes
        if len(time_data) == 0:
            time_data = data[:, 0]

        real_idx = 1
        for _ in range(2 * l + 1):
            psi4_modes_data.append(data[:, real_idx] + 1j * data[:, real_idx + 1])
            real_idx += 2

    return np.array(time_data), np.array(psi4_modes_data)


def find_file_for_l(l: int) -> str:
    """
    Finds the file with the corresponding ell value in the given directory.

    Args:
        l (int): Ell value to search for.

    Returns:
        str: Path to the found file.

    Raises:
        FileNotFoundError: If no file matching the pattern is found.
    """

    for filename in os.listdir(INPUT_DIR):
        if FILE_PATTERN.replace("[MODE=L]", f"{l}") in filename:
            return os.path.join(INPUT_DIR, filename)
    raise FileNotFoundError(f"File with mode l{l} not found.")


def extract_min_omega_ell2_m2(
    time: NDArray[np.float64], data_m2_l2: NDArray[np.complex128]
) -> float:
    """
    Extracts and saves the phase, amplitude, and omega data for l=m=2 mode from psi4 wave.
    Also fits a quadratic to omega and finds its minimum.

    :param time_arr: Array of time data.
    :param mode_data: Dictionary containing the mode data.
    :return: A tuple with parameters from the fit quadratic to omega (minimum value, vertex, curvature).
    """

    collection = psi4_phase_and_amplitude(time, data_m2_l2)
    angular_frequency = collection[3]
    print("hello")
    if OUTPUT_DIR != "":
        print("world")
        labels = ["Time", "Amplitude", "Cum_Phase", "Angular Freq"]
        filename = f"Rpsi4_r{EXT_RAD:06.1f}_ell2_m2_phase_amp_omega.txt"
        arrays_to_txt(labels, collection, filename, OUTPUT_DIR)
        print(f"Time, Amplitude, Phase, and Omega from l=m=2 data saved to {filename}")
        # plt.plot()
        # plt.show()
    return intercept_of_quadratic_fit(time, angular_frequency)


def psi4_phase_and_amplitude(
    time: NDArray[np.float64], cmplx: NDArray[np.complex128]
) -> Tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    """
    Calculates the amplitude and cumulative phase of a gravitational wave signal.

    :param time: A numpy array containing time values.
    :param cmplx: A numpy array containing the a complex signal.

    :return: A tuple containing four numpy arrays (`time`, `amplitude`, `cumulative_phase`, `cumulative_phase_derivative`).

    :raises ValueError: If the lengths of time, real, and imag arrays are not equal.
    """
    if not len(time) == len(cmplx):
        raise ValueError(
            f"The lengths of time {len(time)} and data {len(cmplx)} arrays must be equal."
        )

    # Calculate the amplitude of the gravitational wave signal.
    phases = np.angle(cmplx)
    cycles = 0
    cum_phase = np.zeros_like(time)
    last_phase = phases[0]

    for ph_idx, ph in enumerate(phases):
        # identify phase wrapping
        if np.abs(ph - last_phase) >= np.pi:
            if ph > 0:
                cycles -= 1
            else:
                cycles += 1
        cum_phase[ph_idx] = ph + 2 * np.pi * cycles
        last_phase = ph

    cum_phase_dt = first_time_derivative(time, cum_phase)
    amplitudes = np.abs(cmplx)
    return time, amplitudes, cum_phase, cum_phase_dt


def intercept_of_quadratic_fit(
    time: NDArray[np.float64], data: NDArray[np.float64]
) -> float:
    """
    Fits a quadratic curve to the data within a specified time range and outputs the intercept value.
    This function is intended for l=m=2 angular frequency data.

    :param interval_start: A float specifying the begining of the sample interval.
    :param time: A numpy array containing time values.
    :param data: A numpy array containing data values corresponding to the time values.

    :return: The absolute value of the quadratic curve evaluated at t=0.

    :raises ValueError: If the lengths of time and data arrays are not equal.
    """
    if len(time) != len(data):
        raise ValueError("The lengths of time and data arrays must be equal.")

    # Re-index, keeping only the intersection between numpy arrays
    time_filtered = time[(EXT_RAD <= time) & (time <= EXT_RAD + INTERVAL)]
    data_filtered = data[(EXT_RAD <= time) & (time <= EXT_RAD + INTERVAL)]

    # Fit a quadratic curve to the data using nonlinear least squares
    params, *_ = curve_fit(quadratic, time_filtered, data_filtered)

    # Find the extremum value of the quadratic curve
    a, b, c = params
    # extremum_x = -b / (2 * a)
    # quad_fit_extremum = float(quadratic(extremum_x, a, b, c))
    # print(
    #     f"Quadratic Vertex at (time = {extremum_x:.7e}, value = {quad_fit_extremum:.7e}).\n"
    #     f"Params: a = {a:.7e}, b = {b:.7e}, c = {c:.7e}, Intercept magnitude: {np.fabs(c):.7e}"
    # )
    return np.fabs(c)


def first_time_derivative(
    time: np.ndarray[np.float64],
    data: np.ndarray[np.float64],
) -> np.ndarray[np.float64]:
    """
    Calculate the time derivative of the input data using a second-order finite difference stencil.

    Args:
        time: A numpy array containing time values.
        data: A numpy array containing the data to be differentiated.

    Returns:
        A numpy array containing the time derivative of the input data.

    Raises:
        ValueError: If the lengths of the time and data arrays are not equal.

    >>> time = np.array([0, 1, 2, 3, 4], dtype=np.float64)
    >>> data = np.array([0, 1, 4, 9, 16], dtype=np.float64)
    >>> first_time_derivative(time, data)
    array([1., 2., 4., 6., 7.])
    """

    if len(time) != len(data):
        raise ValueError("Time and data arrays must have the same length.")

    delta_t = time[1] - time[0]
    data_dt = np.zeros_like(data)

    # Second-order in the interior
    data_dt[1:-1] = (data[2:] - data[:-2]) / (2 * delta_t)

    # Drop to first-order at the endpoints
    data_dt[0] = (data[1] - data[0]) / delta_t
    data_dt[-1] = (data[-1] - data[-2]) / delta_t

    return data_dt


def second_time_derivative(
    time: NDArray[np.float64], data: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Compute the second time derivative of the input data using the second-order finite difference method,
    with upwind/downwind stencils for the endpoints.

    :param time: A numpy array containing time values.
    :param data: A numpy array containing data for which the second time derivative is to be calculated.
    :return: A numpy array containing the second time derivative of the input data.

    >>> time = np.array([0, 1, 2, 3, 4], dtype=np.float64)
    >>> data = np.array([0, 1, 4, 9, 16], dtype=np.float64)
    >>> second_time_derivative(time, data)
    array([2., 2., 2., 2., 2.])
    """
    delta_t = time[1] - time[0]
    data_dtdt = np.zeros_like(data)

    # Interior points using central finite difference
    central = data[:-2] - 2 * data[1:-1] + data[2:]
    data_dtdt[1:-1] = central / delta_t**2

    # Endpoint 0: forward finite difference (downwind)
    forward = 2 * data[0] - 5 * data[1] + 4 * data[2] - data[3]
    data_dtdt[0] = forward / delta_t**2

    # Endpoint n-1: backward finite difference (upwind)
    backward = 2 * data[-1] - 5 * data[-2] + 4 * data[-3] - data[-4]
    data_dtdt[-1] = backward / delta_t**2

    return data_dtdt


def quadratic(x: float, a: float, b: float, c: float) -> float:
    """
    Represents a quadratic function.

    :param x: The independent variable.
    :param a: The coefficient of the x^2 term.
    :param b: The coefficient of the x term.
    :param c: The constant term.

    :return: The value of the quadratic function at x.
    """
    return a * x**2 + b * x + c


def get_mode_index(ell: int, em: int, ell_min=ELL_MIN) -> int:
    """
    Returns the array index for mode data given (ell, em).
    The index begins with 0 and increases a first m then l increase.

    Args:
        ell: The l Spherical harmonics mode number
        em: The m Spherical harmonics mode number
        ell_min: The minimum value of ell allowed (default is ELL_MIN).

    Returns:
        The array index for the mode data (ell, em).

    >>> get_mode_index(3, 1)
    9

    """
    return ell**2 + ell + em - ell_min**2


def arrays_to_txt(
    labels: List[str],
    collection: List[NDArray],
    filename: str,
    dir_path: str,
) -> None:
    """Writes an array of NumPy arrays to a text file, formatting each row with labels.

    Args:
        labels: A list of labels, one for each column in the collection.
        collection: A list of NumPy arrays, where each inner array represents a column.
        filename: The name of the file to write to.
        dir_path: The path to the directory where the file will be saved.

    Raises:
        IOError: If there is an error creating the directory or writing to the file.
    """

    try:
        os.makedirs(dir_path, exist_ok=True)
        file_path = os.path.join(dir_path, filename)

        with open(file_path, mode="w", encoding="utf-8") as file:
            file.write(f"# {labels[0]:<12}, ")
            file.write(", ".join([f"{label:<13}" for label in labels[1:]]) + "\n")
            for row in zip(*collection):
                file.write(", ".join([f"{item:.7e}" for item in row]) + "\n")
    except IOError as e:
        raise IOError(f"Error saving data to file: {e}") from e

psi4_fft_to_strain()
