"""
This script generates the gravitational wave meshes to be visualized in VisIt.
It takes all the integrated strain files (suffixed conv_to_strain.txt) for all
the modes and casts them to 3D using linear interpolation and 
spin-weighted spherical harmonic factors. The modes are superimpoed
in the final mesh to recreate the full waveform.
"""

import os, vtk
import time
import numpy as np
import quaternionic, spherical

def read_strain_files(file_path):
    """
    Read an ASCII file with a header describing the real and imaginary parts of the data
    for a specific l mode. Returns the time states and mode data in an easily retrievable format.

    :param file_path: Path to the file to be read.
    :return: A tuple containing the time (numpy array) and a dictionary with keys (l, m) containingth the data.
    :raises ValueError: if the length of the time data is inconsistent across different ell values.  
    """
    mode_data = {}
    time_data_size = -1
    for ell in range(2, 9):
        file_name = file_path.replace("[ELLVAL]", str(ell))
        with open(file_name, "r") as file:
            # Read lines that don't start with '#'
            print(file_name)
            lines = [line for line in file.readlines() if not line.startswith("#")]
    
        # Convert lines to arrays and sort by time
        data = np.array(
            [list(map(np.float64, line.split())) for line in lines]
        )
        data = data[np.argsort(data[:, 0])]

        # Remove duplicate times
        _, index = np.unique(data[:, 0], return_index=True)
        data = data[index]

        # Store time data and check for discrepancies between modes
        time_data = data[:, 0]
        if time_data_size < 0:
            time_data_size = len(time_data)
        elif time_data_size != len(time_data):
            raise ValueError(
                f"Inconsistent time data size for ell={ell}. Expected {time_data_size}, got {len(time_data)}."
            )
        
        # Loop through columns and store real and imaginary parts in dictionary
        for em in range(-ell, ell + 1):
            index = 1 + 2 * (em + ell) # the index for the real valued strain
            mode_data[(ell, em)] = data[:, index] + 1j * data[:, index + 1]

    return time_data, mode_data

def find_swsh_factor(colat, azi, ell, em):
    """
    Calculates the complex valued spin-weighted spherical harmonic (SWSH) factors to be
    multiplied with strain. Uses spherical package to find the factor and the inputted
    colatitude, azimuthal, l, and m values, and a physically defined spin = -2 value.
    
    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param ell: wave mode value l
    :param em: wave mode value m
    :return: a complex valued spin-weighted spherical harmonic factor 
    """
    s = -2
    R = quaternionic.array.from_spherical_coordinates(colat, azi)
    winger = spherical.Wigner(8) # Create a Winger matrix for all modes from l=2 to l=8
    Y = winger.sYlm(s, R)
    swsh_factor = Y[winger.Yindex(ell, em)]
    return swsh_factor

def superimpose_modes_from_angle(colat, azi, num_time_states, mode_data):
    """
    Adds up all the strain modes after factoring in corresponding spin-weighted spherical harmonic
    to specified angle in the mesh. Stored as a new array corresponding to time states.

    :param colat: colatitude angle for the SWSH factor
    :param azi: azimuthal angle for the SWSH factor
    :param num_time_states: number of time states in the strain data
    :param mode_data: dictionary containing strain data for all the modes
    :return: a complex valued numpy array of the superimposed wave
    """
    summation = np.zeros(num_time_states, dtype='complex128')
    for ell in range(2, 9):
        for em in range(-ell, ell + 1):
            swsh_factor = find_swsh_factor(colat, azi, ell, em)
            factored_strain = mode_data[(ell, em)] * swsh_factor
            summation += factored_strain
    return summation

def generate_interpolation_points(time_array, radius_values, r_ext):
    """
    Fills out a 2D array of adjusted time values for the wave strain to be
    linearly interpolated to. First index of the result represents the simulation
    time state (aka which mesh), and the second index represents radial distance to
    interpolate to.

    :param time_array: numpy array of of strain time states.
    :param radius_values: numpy array of the radial points on the mesh.
    :param r_ext: extraction radius of the original data.
    :return: a 2D numpy array of time values.
    """

    time_0 = np.min(time_array)
    time_f = np.max(time_array)
    target_times = np.zeros((len(time_array), len(radius_values)))
    for state, current_time in enumerate(time_array):
        for j, radius in enumerate(radius_values):
            target_time = current_time - radius + r_ext
            if target_time < time_0:
                target_time = time_0
            elif target_time > time_f:
                target_time = time_f
            target_times[state][j] = target_time
    return target_times
    
def initialize_vtk_grid(num_azi, num_radius):
    """
    Sets initial parameters for the mesh generation module and returns
    mesh manipulation objects to write and save data.

    :param num_azi: number of azimuthal points on the mesh
    :param num_radius: number of radial points on the mesh
    :returns: vtk.vtkFloatArray(),
              vtk.vtkUnstructuredGrid(), 
              vtk.vtkPoints(), 
              vtk.vtkXMLUnstructuredGridWriter()
    """
    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    strain_array = vtk.vtkFloatArray()
    strain_array.SetName("Strain")
    strain_array.SetNumberOfComponents(1)
    strain_array.SetNumberOfTuples(num_azi * num_radius)
    cellArray = vtk.vtkCellArray()
    for j in range(num_radius - 1):
        for i in range(num_azi):
            cell = vtk.vtkQuad()
            cell.GetPointIds().SetId(0, i + j * num_azi)
            cell.GetPointIds().SetId(1, (i + 1) % num_azi + j * num_azi)
            cell.GetPointIds().SetId(2, (i + 1) % num_azi + (j + 1) * num_azi)
            cell.GetPointIds().SetId(3, i + (j + 1) * num_azi)
            cellArray.InsertNextCell(cell)
    grid.SetCells(vtk.VTK_QUAD, cellArray)
    writer = (vtk.vtkXMLUnstructuredGridWriter())
    return strain_array, grid, points, writer
    

#######################################

def main():
    """
    Main function that reads the strain data, calculates and factors in spin-weighted spherical harmonics,
    linearly interpolates the strain to fit the mesh points, and creates .vtu mesh file for each time state
    of the simulation. The meshes represent the full superimposed waveform at the polar angle pi/2,
    aka the same plane as the binary black hole merger.
    """
    generic_file_path = r"C:\Users\sethw\OneDrive\PythonFiles\BHVisResearch\r100\Rpsi4_r0100.0_l[ELLVAL]_conv_to_strain.txt"
    output_directory = "./new/meshstates"
    status_messages = True

    # check output directory
    if not os.path.exists(output_directory):
        answer = input(f"Directory: <{output_directory}> does not exist. Would you like to create it (Y/N) ")
        if answer.capitalize() == "Y":
            os.makedirs(output_directory)
        else:
            print("Exiting Program.")
            exit()
    else:
        if os.listdir(output_directory):
            answer = input(f"Data already exists at {output_directory}. Overwrite it? You cannot undo this action. (Y/N) ")
            if answer.capitalize() == "Y":
                for file in os.listdir(output_directory):
                    os.remove(os.path.join(output_directory, file))
            else:
                print("Exiting Program. Change output directory to an empty directory.")
                exit()

    time_array, mode_data = read_strain_files(generic_file_path)
    print(mode_data.keys())
    num_time_states = len(time_array)

    # Mesh generation parameters
    num_radius_points = 450
    num_azi_points = 180
    display_radius = 300
    R_extraction = 100
    amplitude_scale_factor = 600
    omitted_radius_length = 4
    colat = np.pi/2 # colatitude angle representative of the plane of merger

    # Pre-compute theta and radius values for the mesh
    radius_values = np.linspace(0, display_radius, num_radius_points)
    azimuth_values = np.linspace(0, 2 * np.pi, num_azi_points, endpoint=False)
    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)
    
    strain_to_mesh = {} # Holds the final strain points indexed [azimuthal point (key)][time state][radius]
    
    # Apply spin-weighted spherical harmonics, superimpose modes, and interpolate to mesh points
    interpolation_times = generate_interpolation_points(time_array, radius_values, R_extraction)
    for azi_index, azi in enumerate(azimuth_values):
        superimposed_strain = superimpose_modes_from_angle(colat, azi, num_time_states, mode_data)
        strain_to_mesh[azi_index] = np.interp(interpolation_times, time_array, superimposed_strain)
    

    # -----Main Loop: Iterate over all points and construct mesh
    
    strain_array, grid, points, writer = initialize_vtk_grid(num_azi_points, num_radius_points)
    start_time = time.time()
    percentage = np.round(np.linspace(0, num_time_states, 101)).astype(int)

    for state, current_time in enumerate(time_array):
        if status_messages and state == 10:
            end_time = time.time()
            eta = (end_time - start_time) * num_time_states / 10
            print(f"Creating {num_time_states} meshes and saving them to: {output_directory}\nEstimated time: {int(eta / 60)} minutes")
        if status_messages and state != 0 and np.isin(state, percentage):
            print(f" {int(state * 100 / (num_time_states - 1))}% done", end="\r")

        # Define the points and their coordinates
        points.Reset()
        index = 0
        for j, radius in enumerate(radius_values):
            for i, azi in enumerate(azimuth_values):
                h_tR = strain_to_mesh[i][state][j]
                x = x_values[j, i]
                y = y_values[j, i]
                if radius <= omitted_radius_length: # Introduce a discontinuity to make room for the Black Holes
                    z = np.nan
                    strain_value = np.nan
                else:
                    strain_value = h_tR.real
                    z = strain_value * amplitude_scale_factor
                points.InsertNextPoint(x, y, z)
                strain_array.SetTuple1(index, strain_value)
                index += 1

        grid.SetPoints(points)
        grid.GetPointData().AddArray(strain_array)
        writer.SetFileName(output_directory + f"/state{state}.vtu")
        writer.SetInputData(grid)
        writer.Write()

if __name__ == "__main__":
    main()