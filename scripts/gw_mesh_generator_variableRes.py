import vtk
import os
import time
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm

# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test7_variableRes/"

# Set up grid parameters
status_messages = True
plot_strain = True
l = 2
m = 2
x_lim = 200
y_lim = 200
z_lim = 1

# Mesh parameters
resolution = 100
display_radius = 100
R_ext = 100

# Calculate spin-weighted spherical harmonics for every point in the mesh, returns them as complex values in a 2d array
def set_sph_harm_array(l, m, s):
    global resolution, display_radius

    sph_harm_points = np.zeros((resolution, resolution), dtype=np.complex128)
    interval = 2 * display_radius / resolution

    for j in range(resolution):
        y = -(display_radius - 1) + (j * interval)
        for i in range(resolution):
            x = -(display_radius - 1) + (i * interval)
            r = np.sqrt(x ** 2 + y ** 2)
            theta = np.pi / 2
            phi = np.arctan2(y, x)
            Y_lm = sph_harm(l, m, phi, theta)
            Y = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * np.math.factorial(l - m) / np.math.factorial(l + m)) * Y_lm
            sph_harm_points[i, j] = Y

    return sph_harm_points

# This function allows a linearly interpolated strain to be found given a target time
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

    sph_harm_points = set_sph_harm_array(2, 2, -2)

    return length, h_strain, sph_harm_points, h_time

length, h_strain, sph_harm_points, h_time = initialize()
def show_strain_plot():
    plt.plot(h_time, h_strain.real)
    plt.title("Gravitational Wave Strain vs Time\nl =" + str(l) + "m =" + str(m) + "\nData Extraction Radius =" + str(
        R_ext) + "meters")
    plt.xlabel("Time")
    plt.ylabel("Re Strain")

if plot_strain:
    show_strain_plot()

start_time = time.time()
percentage = np.round(np.linspace(0, length, 101)).astype(int)

state = 0
for current_time in h_time:
    state += 1
    points = vtk.vtkPoints()
    grid = vtk.vtkStructuredGrid()
    
    t = np.where(h_time == current_time)[0][0]
    if status_messages and t == 10:
        end_time = time.time()
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}")
    if status_messages and t != 0 and np.isin(t, percentage):
        print(f" {int(t * 100 / (length - 1))}% done", end="\r")

    grid.SetDimensions(resolution, resolution, z_lim)

    highRes_factor = 6
    lowRes_factor = 2
    highRes_bound = 20
    for current_dr in range(display_radius):
        if current_dr <= highRes_bound:
            resolution = highRes_bound * highRes_factor 
            interval = 2 * display_radius / resolution
            
            def set_sph_harm_array(l, m, s):
                global resolution, display_radius

                sph_harm_points = np.zeros((resolution, resolution), dtype=np.complex128)
                interval = 2 * display_radius / resolution

                for j in range(resolution):
                    y = -(display_radius - 1) + (j * interval)
                    for i in range(resolution):
                        x = -(display_radius - 1) + (i * interval)
                        r = np.sqrt(x ** 2 + y ** 2)
                        theta = np.pi / 2
                        phi = np.arctan2(y, x)
                        Y_lm = sph_harm(l, m, phi, theta)
                        Y = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * np.math.factorial(l - m) / np.math.factorial(l + m)) * Y_lm
                        sph_harm_points[i, j] = Y
                return sph_harm_points
            sph_harm_points = set_sph_harm_array(2, 2, -2)

            for j in range(resolution):
                y = -(display_radius - 1) + (j * interval)
                for i in range(resolution):
                    x = -(display_radius - 1) + (i * interval)
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

                    scale_factor = 600
                    z = (Y.real * h_tR.real - Y.imag * h_tR.imag) * scale_factor
                    points.InsertNextPoint(x, y, z)
        else:
            for j in range(resolution):
                y = -(display_radius - 1) + (j * interval)
                for i in range(resolution):
                    x = -(display_radius - 1) + (i * interval)
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

                    scale_factor = 600
                    z = (Y.real * h_tR.real - Y.imag * h_tR.imag) * scale_factor
                    points.InsertNextPoint(x, y, z)

    grid.SetPoints(points)

    # Create strain data array
    strain_array = vtk.vtkFloatArray()
    strain_array.SetName("Strain")
    strain_array.SetNumberOfComponents(1)
    strain_array.SetNumberOfTuples(resolution * resolution)

    # Populate strain data array
    index = 0
    for j in range(resolution):
        for i in range(resolution):
            x = -(display_radius - 1) + (i * interval)
            y = -(display_radius - 1) + (j * interval)
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

            strain_value = (Y.real * h_tR.real - Y.imag * h_tR.imag) * scale_factor
            strain_array.SetTuple1(index, strain_value)
            index += 1

    grid.GetPointData().AddArray(strain_array)

    writer = vtk.vtkXMLStructuredGridWriter()
    filename = output_directory + f"/state{state}.vts"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()

print("Mesh database completed in", output_directory)