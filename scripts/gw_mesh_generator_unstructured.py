import vtk
import os
import math
import time
import numpy as np
from scipy.special import sph_harm

# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test8_polar_zeroR/"

# Define the dimensions of the grid
numTheta = 180  # Number of points along the theta direction
numRadius = 450 # Number of points along the radius direction
new_numTheta = 720
new_numRadius = 300
change_time = 4150
display_radius = 200 # Mesh radius
R_ext = 100 # Extraction radius

# Calculate spin-weight spherical harmonic
def set_sph_harm_array(l, m, s, numTheta, numRadius, display_radius):
    # Initialize the spherical harmonic array
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
            Y = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * np.math.factorial(l - m) / np.math.factorial(l + m)) * Y_lm

            # Store the spherical harmonic in the array
            sph_harm_points[i, j] = Y

    return sph_harm_points



# Linearly interpolate strain for any given time
def interpolated_strain(target_time, source_time, data):
    if len(source_time) != len(data):
        raise ValueError("source_time and data must have the same number of rows")
    if not (np.diff(source_time) > 0).all():
        raise ValueError("source_time must be strictly increasing")

    interpolated_data = np.interp(target_time, source_time, data)
    return interpolated_data

# Reads inputted strain data - returns data file row length, initial strain data (complex), spin-weighted spherical harmonics (complex), and time values
def initialize():
    if os.path.exists(output_directory):
        if len(os.listdir(output_directory)) != 0:
            answer = input(f"Data already exists at {output_directory}. Overwrite it? You cannot undo this action. (Y/N) ")
            if answer.capitalize() == "Y":
                for file in os.listdir(output_directory):
                    os.remove(f"{output_directory}/{file}")
            else:
                print("Exiting Program. Change output directory to an empty directory.")
                exit()
    else:
        last_slash_index = output_directory.rfind("/")
        super_directory = output_directory[:last_slash_index] if last_slash_index != -1 else output_directory
        if os.path.exists(super_directory):
            os.makedirs(output_directory)
        else:
            print(f"Error: {super_directory} does not exist.")
            exit()

    def valid_line(line):
        return not line.startswith("#")

    with open(input_file, 'r') as f:
        strain_data = np.array([list(map(float, line.split())) for line in f if valid_line(line)])

    strain_data = np.unique(strain_data, axis=0)
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j * h_imag

    length = len(h_strain)

    sph_harm_points = set_sph_harm_array(2, 2, -2, numTheta, numRadius, display_radius)

    return length, h_strain, sph_harm_points, h_time

length, h_strain, sph_harm_points, h_time = initialize()

# Iterate over all points and construct mesh
start_time = time.time()
percentage = np.round(np.linspace(0, length, 101)).astype(int)
state = 0
scale_factor = 600
status_messages = True
omit_radius = 3

# Main Loop
for current_time in h_time:
    state += 1
    t = np.where(h_time == current_time)[0][0]
    if status_messages and t == 10:
        end_time = time.time()
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}")
    if status_messages and t != 0 and np.isin(t, percentage):
        print(f" {int(t * 100 / (length - 1))}% done", end="\r")
    '''
    # Change resolution at the specified time
    if current_time >= change_time:
        numTheta = new_numTheta
        numRadius = new_numRadius
        sph_harm_points = set_sph_harm_array(2, 2, -2, numTheta, numRadius, display_radius)
    '''
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
    for j in range(numRadius):
        for i in range(numTheta):
            theta = 2.0 * vtk.vtkMath.Pi() * float(i) / float(numTheta)
            radius = display_radius * float(j) / float(numRadius - 1)  # Scaling radius from 0 to 1
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            current_r = np.sqrt(x ** 2 + y ** 2)
            target_time = current_time - current_r + R_ext
            time_0 = np.min(h_time)
            time_f = np.max(h_time)
            if target_time < time_0:
                target_time = time_0
            elif target_time > time_f:
                target_time = time_f
            h_tR = interpolated_strain(target_time, h_time, h_strain)
            Y = sph_harm_points[i, j]
            z = (Y.real * h_tR.real - Y.imag * h_tR.imag) * scale_factor
            strain_value = z
            if current_r <= 4:
                z = -10
                strain_value = 0
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

    writer = vtk.vtkXMLUnstructuredGridWriter() #vtk.vtkXMLPUnstructuredGridWriter() produces paralell files
    filename = output_directory + f"/state{state}.vtu"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()