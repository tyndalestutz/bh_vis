import os, math, time, vtk
import numpy as np
from scipy.special import sph_harm

# File parameters
input_file = "./Rpsi4_l2-r0100.0.txt"
output_directory = "./outputmeshstates/"
status_messages = True

# Plot parameters
numRadius = 450  # Number of points along the radius direction
numTheta = 180  # Number of points along the theta direction
display_radius = 300  # Mesh radius
R_ext = 100  # Extraction radius
l = 2
m = 2
s = -2
scale_factor = 600
omit_radius = 4
plot_strain = True

def set_sph_harm_array(l, m, s, radius_array, theta_array): #uses scipy.special sph_harm
    sph_harm_points = np.zeros((len(radius_array), len(theta_array)), dtype=np.complex128)
    phi = np.pi / 2  # phi is fixed to pi/2, why?
    spin_weight = (-1) ** s * np.sqrt((2 * l + 1) / (4 * np.pi) * math.factorial(l - m) / math.factorial(l + m)) # why?
    for j, radius in enumerate(radius_array):
        for i, theta in enumerate(theta_array):
            sph_harm_points[j, i] = spin_weight * sph_harm(m, l, theta, phi)
    return sph_harm_points

# Reads inputted strain data - returns data file row length, initial strain data (complex), and time values
def initialize():
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

    with open(input_file, "r") as f:
        strain_data = np.array([list(map(float, line.split())) for line in f if not(line.startswith("#"))])

    strain_data = np.unique(strain_data, axis=0)
    length = len(strain_data)
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j * h_imag
    if not (np.diff(h_time) > 0).all():
        raise ValueError("h_time must be strictly increasing")
    return length, h_strain, h_time


length, h_strain, h_time = initialize()
percentage = np.round(np.linspace(0, length, 101)).astype(int)

# Pre-compute the theta and radius values
radius_values = np.linspace(0, display_radius, numRadius)
theta_values = np.linspace(0, 2 * np.pi, numTheta, endpoint=False)
rv, tv = np.meshgrid(radius_values, theta_values, indexing="ij")
x_values = rv * np.cos(tv)
y_values = rv * np.sin(tv)
sph_array = set_sph_harm_array(l, m, s, radius_values, theta_values)

#generate and filter target times to pre-interpolate
time_0 = np.min(h_time)
time_f = np.max(h_time)
t_array = np.zeros((length, numRadius))
for state, current_time in enumerate(h_time):
    for j, radius in enumerate(radius_values):
        target_time = current_time - radius + R_ext
        if target_time < time_0:
            target_time = time_0
        elif target_time > time_f:
            target_time = time_f
        t_array[state][j] = target_time
h_array = np.interp(t_array, h_time, h_strain)

# Intitialize vtk grid, points, data-array, grid cells, and writer
grid = vtk.vtkUnstructuredGrid()
points = vtk.vtkPoints()
strain_array = vtk.vtkFloatArray()
strain_array.SetName("Strain")
strain_array.SetNumberOfComponents(1)
strain_array.SetNumberOfTuples(numTheta * numRadius)
cellArray = vtk.vtkCellArray()
for j in range(numRadius - 1):
    for i in range(numTheta):
        cell = vtk.vtkQuad()
        cell.GetPointIds().SetId(0, i + j * numTheta)
        cell.GetPointIds().SetId(1, (i + 1) % numTheta + j * numTheta)
        cell.GetPointIds().SetId(2, (i + 1) % numTheta + (j + 1) * numTheta)
        cell.GetPointIds().SetId(3, i + (j + 1) * numTheta)
        cellArray.InsertNextCell(cell)
grid.SetCells(vtk.VTK_QUAD, cellArray)
writer = (vtk.vtkXMLUnstructuredGridWriter())

# Main Loop: Iterate over all points and construct mesh
start_time = time.time()
for state, current_time in enumerate(h_time):
    if status_messages and state == 10:
        end_time = time.time()
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to: {output_directory}\nEstimated time: {int(eta / 60)} minutes")
    if status_messages and state != 0 and np.isin(state, percentage):
        print(f" {int(state * 100 / (length - 1))}% done", end="\r")

    # Define the points and their coordinates
    points.Reset()
    index = 0
    for j, radius in enumerate(radius_values):
        h_tR = h_array[state][j]
        for i, theta in enumerate(theta_values):
            x = x_values[j, i]
            y = y_values[j, i]
            S = sph_array[j, i]
            if radius <= omit_radius: # Introduce a discontinuity to make room for the Black Holes
                z = np.nan
                strain_value = np.nan
            else:
                strain_value = S.real * h_tR.real - S.imag * h_tR.imag
                z = strain_value * scale_factor
            points.InsertNextPoint(x, y, z)
            strain_array.SetTuple1(index, strain_value)
            index += 1

    grid.SetPoints(points)
    grid.GetPointData().AddArray(strain_array)

    writer.SetFileName(output_directory + f"/state{state}.vtu")
    writer.SetInputData(grid)
    writer.Write()