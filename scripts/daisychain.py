import sys
import os
import time as mod_time
from typing import Dict, List, Tuple, Any

# import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.optimize import curve_fit
import quaternionic
import spherical
import vtk

OUTPUT_DIR = "visit/dataout"
VTKOUT = "visit/output"
INPUT_DIR = "r100"
FILE_PATTERN = "_l[MODE=L]-"
ELL_MIN = 2
ELL_MAX = 8
S_MODE = -2
EXT_RAD = 100
INTERVAL = 200
CUTTOFF_FACTOR = 0.75
STATUS_MESSAGES = True

def find_file_for_l(ell: int) -> str:
    """
    Finds the file path with the corresponding ell value in the given directory.

    :param ell: (int): l mode to search for.
    :return: Path to the found file.
    :raises FileNotFoundError: If no file matching the pattern is found.
    """

    for filename in os.listdir(INPUT_DIR):
        if FILE_PATTERN.replace("[MODE=L]", f"{ell}") in filename:
            return os.path.join(INPUT_DIR, filename)
    raise FileNotFoundError(f"File with mode ell{ell} not found.")


def read_psi4_dir() -> tuple[np.ndarray, np.ndarray]:
    """
    Reads data from psi4 output directory and returns time and mode data.

    :return: tuple[np.ndarray, np.ndarray]
        - time_data: Array of numpy.float64 time values (shape: (n_times,) ).
        - mode_data: 2D Array for modes of numpy.complex128 data (shape: (2*l+1, n_times,) ).
    """

    time_data: list[np.ndarray] = []
    psi4_modes_data: list[np.ndarray] = []
    consistant_length = -1

    for l in range(ELL_MIN, ELL_MAX + 1):
        filepath = find_file_for_l(l)
        with open(filepath, "r", encoding="utf-8") as file:
            lines = [line for line in file.readlines() if not line.startswith("#")]
        data = np.array([np.array(line.split(), dtype=np.float64) for line in lines])

        time_data, index = np.unique(
            data[:, 0], return_index=True
        )  # sorts by time, removing duplicates
        data = data[index]  # sorts data accordningly

        if consistant_length == -1:
            consistant_length = len(time_data)

        if consistant_length != len(time_data):
            raise ValueError(
                f"Inconsistent time data for ell={l}. Expected {consistant_length}, got {len(time_data)}."
            )

        real_idx = 1
        for _ in range(2 * l + 1):
            psi4_modes_data.append(data[:, real_idx] + 1j * data[:, real_idx + 1])
            real_idx += 2
    return np.array(time_data), np.array(psi4_modes_data)


def get_index_from_modes(ell: int, em: int, ell_min=ELL_MIN) -> int:
    """
    Returns the array index for mode data given (ell, em).
    The index begins with 0 and through m (inner loop) then l (outer loop).

    :param ell: The l Spherical harmonics mode number
    :param em: The m Spherical harmonics mode number
    :param ell_min: The minimum ell value used in the array (default is ELL_MIN).

    :return: The mode data array index for (ell, em).

    >>> get_index_from_modes(3, 1, 2)
    9

    """
    return ell**2 + ell + em - ell_min**2


def get_modes_from_index(idx: int, ell_min=ELL_MIN) -> int:
    idx += ell_min**2
    ell = int(np.sqrt(idx))
    em = idx - ell**2 - ell
    return ell, em


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


def first_time_derivative(
    time: NDArray[np.float64],
    data: NDArray[np.float64],
) -> NDArray[np.float64]:
    """
    Calculate the time derivative of the input data using a second-order finite difference stencil.

    :param time: A numpy array containing time values.
    :param data: A numpy array containing the data to be differentiated.
    :return: A numpy array containing the time derivative of the input data.

    >>> time = np.array([0, 1, 2, 3, 4], dtype=np.float64)
    >>> data = np.array([0, 1, 4, 9, 16], dtype=np.float64)
    >>> first_time_derivative(time, data)
    array([1., 2., 4., 6., 7.])
    """

    delta_t = time[1] - time[0]
    data_dt = np.zeros_like(data)

    # Second-order in the interior:
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


def psi4_phase_and_amplitude(
    time: NDArray[np.float64], cmplx: NDArray[np.complex128]
) -> Tuple[
    NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]
]:
    """
    Calculates the amplitude and cumulative phase of a gravitational wave signal.

    :param time: A numpy array containing time values.
    :param cmplx: A numpy array containing the a complex signal.
    :return: A tuple containing four numpy arrays (time, amplitude, cumulative_phase, cumulative_phase_derivative).
    :raises ValueError: If the lengths of time, real, and imag arrays are not equal.
    """
    if not len(time) == len(cmplx):
        raise ValueError(
            f"The lengths of time {len(time)} and data {len(cmplx)} arrays must be equal."
        )

    amplitudes = np.abs(cmplx)
    phases = np.angle(cmplx)
    cycles = 0
    cum_phase = np.zeros_like(time)
    last_phase = phases[0]

    for i, ph in enumerate(phases):
        # identify phase wrapping
        if np.abs(ph - last_phase) >= np.pi:
            cycles += -1 if ph > 0 else 1
        cum_phase[i] = ph + 2 * np.pi * cycles
        last_phase = ph

    cum_phase_dt = first_time_derivative(time, cum_phase)
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

    def quadratic(x: float, a: float, b: float, c: float) -> float:
        """
        Evaluates a quadratic (ax^2 + bx + c)

        :param x: The independent variable.
        :param a: The coefficient of the x^2 term.
        :param b: The coefficient of the x term.
        :param c: The constant term.

        :return: The value of the quadratic function at x.
        """
        return a * x**2 + b * x + c

    # Re-index, keeping only the intersection between numpy arrays
    time_filtered = time[(EXT_RAD <= time) & (time <= EXT_RAD + INTERVAL)]
    data_filtered = data[(EXT_RAD <= time) & (time <= EXT_RAD + INTERVAL)]

    # Fit a quadratic curve to the data using nonlinear least squares
    params, *_ = curve_fit(quadratic, time_filtered, data_filtered)

    # Find the extremum value of the quadratic curve
    a, b, c = params
    extremum_x = -b / (2 * a)
    quad_fit_extremum = float(quadratic(extremum_x, a, b, c))
    print(
        f"Quadratic Vertex at (time = {extremum_x:.7e}, value = {quad_fit_extremum:.7e}).\n"
        f"Params: a = {a:.7e}, b = {b:.7e}, c = {c:.7e}, Intercept magnitude: {np.fabs(c):.7e}"
    )
    return np.fabs(c)


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
    if OUTPUT_DIR != "":
        labels = ["Time", "Amplitude", "Cum_Phase", "Angular Freq"]
        filename = f"Rpsi4_r{EXT_RAD:06.1f}_ell2_m2_phase_amp_omega.txt"
        arrays_to_txt(labels, collection, filename, OUTPUT_DIR)
        print(f"Time, Amplitude, Phase, and Omega from l=m=2 data saved to {filename}")
        # plt.plot()
        # plt.show()
    return intercept_of_quadratic_fit(time, angular_frequency)


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
        time_arr, psi4_modes_data[get_index_from_modes(2, 2)]
    )
    freq_cutoff = CUTTOFF_FACTOR * min_omega_l2m2

    # Initialize arrays for strain modes and their second time derivatives
    strain_modes = np.zeros_like(psi4_modes_data)
    strain_modes_ddot = np.zeros_like(psi4_modes_data)
    labels = ["time"]

    # Calculate frequency list for FFT
    freq_list = np.fft.fftfreq(len(time_arr), time_arr[1] - time_arr[0]) * 2 * np.pi

    # Next loop over modes and perform an FFT:
    mode_idx = 0
    for l in range(ELL_MIN, ELL_MAX + 1):
        for m in range(-l, l + 1):
            # Apply FFT and filter, see Eq. 27 in https://arxiv.org/abs/1006.1632
            fft_result = np.fft.fft(psi4_modes_data[get_index_from_modes(l, m)])
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
            # plt.plot(time_arr, mode_data[mode_idx] - strain_modes_ddot[mode_idx])
            labels.append(f"(l={l} m={m})")
            mode_idx += 1
        # plt.show()
    # Save the strain output to a file with _conv_to_strain.txt extension
    arrays_to_txt(labels, np.concatenate(([time_arr],strain_modes)), "psi4_conv_to_strain.csv", OUTPUT_DIR)
    print(f"Unweighted Strain Modes saved to {OUTPUT_DIR}/psi4_conv_to_strain.csv")
    return time_arr, strain_modes


def swsh_summation_angles(colat: float, azi: NDArray[np.float64], mode_data):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as an array corresponding to [angle, time] time_idxs.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param n_times: number of time time_idxs in the strain data
    :param mode_data: dictionary containing strain data for all the modes
    :return: a complex valued numpy array of the superimposed wave
    """
    quat_arr = quaternionic.array.from_spherical_coordinates(colat, azi)
    winger = spherical.Wigner(ELL_MAX, ELL_MIN)
    # Create an swsh array shaped like (n_modes, n_quaternions)
    swsh_arr = winger.sYlm(S_MODE, quat_arr).T
    # mode_data has shape (n_modes, n_times), swsh_arr has shape (n_mode, n_pts).
    # Pairwise multiply and sum over modes: the result has shape (n_pts, n_times).
    return np.sum(mode_data[:, np.newaxis, :] * swsh_arr[:, :, np.newaxis], axis=0)


def generate_interpolation_points(  # Could use some revisiting, currently keeps n_times constant
    time_array: NDArray[np.float64], radius_values: NDArray[np.float64], r_ext: int
) -> NDArray[np.float64]:
    """
    Fills out a 2D array of adjusted time values for the wave strain to be
    linearly interpolated to. First index of the result represents the simulation
    time time_idx (aka which mesh), and the second index represents radial distance to
    interpolate to.

    :param time_array: numpy array of of strain time time_idxs.
    :param radius_values: numpy array of the radial points on the mesh.
    :param r_ext: extraction radius of the original data.
    :return: a 2D numpy array (n_radius, n_times) of time values.
    """
    # Repeat time_array and radius_values to match the shape of the output array
    time_repeated = np.repeat(time_array[np.newaxis, :], len(radius_values), axis=0)
    radius_repeated = np.repeat(radius_values[:, np.newaxis], len(time_array), axis=1)

    target_times = (
        time_repeated - radius_repeated + r_ext
    )  # Shape is (n_radius, n_times)
    return np.clip(target_times, time_array.min(), time_array.max())


def initialize_vtk_grid(n_azi: int, n_radius: int):
    """
    Sets initial parameters for the mesh generation module and returns
    mesh manipulation objects to write and save data.

    :param n_azi: number of azimuthal points on the mesh
    :param n_radius: number of radial points on the mesh
    :returns: vtk.vtkFloatArray(),
              vtk.vtkUnstructuredGrid(),
              vtk.vtkPoints(),
              vtk.vtkXMLUnstructuredGridWriter()
    """
    grid = vtk.vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    data_vtk_arr = vtk.vtkFloatArray()
    data_vtk_arr.SetName("Strain")
    data_vtk_arr.SetNumberOfComponents(1)
    data_vtk_arr.SetNumberOfTuples(n_azi * n_radius)
    cell_array = vtk.vtkCellArray()
    for j in range(n_radius - 1):
        for i in range(n_azi):
            cell = vtk.vtkQuad()
            cell.GetPointIds().SetId(0, i + j * n_azi)
            cell.GetPointIds().SetId(1, (i + 1) % n_azi + j * n_azi)
            cell.GetPointIds().SetId(2, (i + 1) % n_azi + (j + 1) * n_azi)
            cell.GetPointIds().SetId(3, i + (j + 1) * n_azi)
            cell_array.InsertNextCell(cell)
    grid.SetCells(vtk.VTK_QUAD, cell_array)
    writer = vtk.vtkXMLUnstructuredGridWriter()
    return data_vtk_arr, grid, pts, writer


def check_output_directory(output_directory):
    """
    Creates the output directory for the .vtu files,
    and asks to overwrite existing files if they exist.

    :param output_directory: directory for function to create.
    :returns: None
    """
    os.makedirs(output_directory, exist_ok=True)
    if os.listdir(output_directory):
        answer = input(
            f"""Data already exists at {output_directory}.\nOverwrite it? You cannot undo this action. (Y/N): """
        )
        if answer.capitalize() == "Y":
            for file in os.listdir(output_directory):
                os.remove(os.path.join(output_directory, file))
        else:
            print("Exiting Program. Change output directory to an empty directory.")
            sys.exit()


def gen_wave_meshes():
    output_directory = VTKOUT
    check_output_directory(output_directory)

    # Mesh generation parameters
    ext_rad = EXT_RAD
    n_rad_pts = 450
    n_azi_pts = 180
    display_radius = 300
    amplitude_scale_factor = 600
    omitted_radius_length = 4
    colat = np.pi / 2  # This colatitude angle represents the BH merger plane

    time_array, mode_data = psi4_fft_to_strain()
    n_times = len(time_array)

    # Pre-compute theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, n_rad_pts)
    azimuth_values = np.linspace(0, 2 * np.pi, n_azi_pts, endpoint=False)

    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)

    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    strain_to_mesh = np.zeros((n_rad_pts, n_azi_pts, n_times))

    strain_azi = swsh_summation_angles(colat, azimuth_values, mode_data).real

    lerp_times = generate_interpolation_points(time_array, radius_values, ext_rad)
    for i in range(n_azi_pts):
        strain_to_mesh[:, i, :] = np.interp(lerp_times, time_array, strain_azi[i, :])

    # -----Main Loop: Iterate over all points and construct mesh

    strain_vtk_arr, vtk_grid, vtk_pts, writer = initialize_vtk_grid(
        n_azi_pts, n_rad_pts
    )

    if STATUS_MESSAGES:
        start_time = mod_time.time()
        percentage = np.round(np.linspace(0, n_times, 101)).astype(int)

    for time_idx in range(n_times):
        if STATUS_MESSAGES:
            if time_idx == 10:
                end_time = mod_time.time()
                eta = (end_time - start_time) * n_times / 10
                print(
                    f"""Creating {n_times} meshes and saving them to: {output_directory}\nEstimated time: {int(eta / 60)} minutes"""
                )

            if time_idx != 0 and np.isin(time_idx, percentage):
                print(f" {int(time_idx * 100 / (n_times - 1))}% done", end="\r")

        # Define the points and their coordinates
        vtk_pts.Reset()
        index = 0
        for rad_idx, radius in enumerate(radius_values):
            for azi_idx in range(n_azi_pts):
                x = x_values[rad_idx, azi_idx]
                y = y_values[rad_idx, azi_idx]
                # Introduces a discontinuity to make room for the Black Holes
                if radius <= omitted_radius_length:
                    strain_value = np.nan
                else:
                    strain_value = strain_to_mesh[rad_idx][azi_idx][time_idx]
                z = strain_value * amplitude_scale_factor
                vtk_pts.InsertNextPoint(x, y, z)
                strain_vtk_arr.SetTuple1(index, strain_value)
                index += 1

        vtk_grid.SetPoints(vtk_pts)
        vtk_grid.GetPointData().AddArray(strain_vtk_arr)
        writer.SetFileName(os.path.join(output_directory, f"state{time_idx}.vtu"))
        writer.SetInputData(vtk_grid)
        writer.Write()


if __name__ == "__main__":
    # First run doctests
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Then run the main() function.
    gen_wave_meshes()