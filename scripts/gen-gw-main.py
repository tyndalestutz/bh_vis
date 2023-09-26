import vtk
import os
import math
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm

# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "./gw_test12/"

# Parameters
numRadius = 450  # Number of points along the radius direction
numTheta = 180  # Number of points along the theta direction
display_radius = 300  # Mesh radius
R_ext = 100  # Extraction radius
l = 2
m = 2
s = -2
plot_strain = True


def set_sph_harm_array(l, m, s, radius_array, theta_array):
    # Initialize the spherical harmonic array
    sph_harm_points = np.zeros(
        (len(radius_array), len(theta_array)), dtype=np.complex128
    )
    phi = np.pi / 2  # phi is fixed to pi/2, similar to the original script
    yspin = (-1) ** s * np.sqrt(
        (2 * l + 1) / (4 * np.pi) * np.math.factorial(l - m) / np.math.factorial(l + m)
    )
    # Loop over all points in the grid
    for j, radius in enumerate(radius_array):
        for i, theta in enumerate(theta_array):
            # Calculate the radial and angular coordinates
            # Calculate the spherical harmonic
            Y = yspin * sph_harm(m, l, theta, phi)
            # Adjust for spin-weighted
            # Store the spherical harmonic in the array
            sph_harm_points[j, i] = Y
    return sph_harm_points


# Linearly interpolate strain for any given time
def interpolated_strain(target_time, source_time, data):
    if len(source_time) != len(data):
        raise ValueError("source_time and data must have the same number of rows")
    if not (np.diff(source_time) > 0).all():
        raise ValueError("source_time must be strictly increasing")

    interpolated_data = np.interp(target_time, source_time, data)
    return interpolated_data


# Reads inputted strain data - returns data file row length, initial strain data (complex), and time values
def initialize():
    if not os.path.exists(output_directory):
        super_directory = os.path.dirname(output_directory)
        if os.path.exists(super_directory):
            os.makedirs(output_directory)
        else:
            print(f"Error: {super_directory} does not exist.")
            exit()
    else:
        if os.listdir(output_directory):
            answer = input(
                f"Data already exists at {output_directory}. Overwrite it? You cannot undo this action. (Y/N) "
            )
            if answer.capitalize() == "Y":
                for file in os.listdir(output_directory):
                    os.remove(os.path.join(output_directory, file))
            else:
                print("Exiting Program. Change output directory to an empty directory.")
                exit()

    def valid_line(line):
        return not line.startswith("#")

    with open(input_file, "r") as f:
        strain_data = np.array(
            [list(map(float, line.split())) for line in f if valid_line(line)]
        )

    strain_data = np.unique(strain_data, axis=0)
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j * h_imag

    length = len(h_strain)

    return length, h_strain, h_time


length, h_strain, h_time = initialize()

# Iterate over all points and construct mesh

percentage = np.round(np.linspace(0, length, 101)).astype(int)
scale_factor = 600
status_messages = True
omit_radius = 4

# Pre-compute the theta and radius values outside the loop
radius_values = np.linspace(0, display_radius, numRadius)
theta_values = np.linspace(0, 2 * np.pi, numTheta, endpoint=False)
rv, tv = np.meshgrid(radius_values, theta_values, indexing="ij")
x_values = rv * np.cos(tv)
y_values = rv * np.sin(tv)
sph_harm_points = set_sph_harm_array(l, m, s, radius_values, theta_values)
time_0 = np.min(h_time)
time_f = np.max(h_time)

t_array = np.zeros((len(h_time), len(radius_values)))
for state, current_time in enumerate(h_time, start=1):
    for j, radius in enumerate(radius_values):
        target_time = current_time - radius + R_ext
        if target_time < time_0:
            target_time = time_0
        elif target_time > time_f:
            target_time = time_f
        t_array[state][j] = target_time
h_array = np.interp(t_array, h_time, h_strain)

start_time = time.time()
# Main Loop
for state, current_time in enumerate(h_time, start=1):
    if status_messages and state == 11:
        end_time = time.time()
        eta = (end_time - start_time) * length / 10
        print(
            f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}"
        )
    if status_messages and state != 0 and np.isin(state, percentage):
        print(f" {int(state * 100 / (length - 1))}% done", end="\r")
    """
    # Change resolution at the specified time
    if current_time >= change_time:
        numTheta = new_numTheta
        numRadius = new_numRadius
        sph_harm_points = set_sph_harm_array(l, m, s, numTheta, numRadius, display_radius)
    """
    # Create vtkUnstructuredGrid
    grid = vtk.vtkUnstructuredGrid()

    # Create strain data array
    strain_array = vtk.vtkFloatArray()
    strain_array.SetName("Strain")
    strain_array.SetNumberOfComponents(1)
    strain_array.SetNumberOfTuples(numTheta * numRadius)

    # Define the points and their coordinates
    points = vtk.vtkPoints()
    index = 0
    for j, radius in enumerate(radius_values):
        h_tR = h_array[state][j]
        for i, theta in enumerate(theta_values):
            x = x_values[j, i]
            y = y_values[j, i]
            Y = sph_harm_points[j, i]
            strain_value = Y.real * h_tR.real - Y.imag * h_tR.imag
            z = strain_value * scale_factor
            # Introduce a discontinuity to make room for the Black Holes
            if radius <= omit_radius:
                z = np.nan
                strain_value = np.nan
            points.InsertNextPoint(x, y, z)
            strain_array.SetTuple1(index, strain_value)
            index += 1

    # Set the points in the vtkUnstructuredGrid and strain_array
    grid.SetPoints(points)
    grid.GetPointData().AddArray(strain_array)

    # Define the connectivity of the cells
    cellArray = vtk.vtkCellArray()
    for j in range(numRadius - 1):
        for i in range(numTheta):
            cell = vtk.vtkQuad()
            cell.GetPointIds().SetId(0, i + j * numTheta)
            cell.GetPointIds().SetId(1, (i + 1) % numTheta + j * numTheta)
            cell.GetPointIds().SetId(2, (i + 1) % numTheta + (j + 1) * numTheta)
            cell.GetPointIds().SetId(3, i + (j + 1) * numTheta)
            cellArray.InsertNextCell(cell)

    # Set the cells in the vtkUnstructuredGrid
    grid.SetCells(vtk.VTK_QUAD, cellArray)

    writer = (
        vtk.vtkXMLUnstructuredGridWriter()
    )  # vtk.vtkXMLPUnstructuredGridWriter() produces paralell files
    filename = output_directory + f"/state{state}.vtu"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()
