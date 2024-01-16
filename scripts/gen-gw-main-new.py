import os, time, vtk, numpy as np
import re

# Modified from gen-gw-main.py to use new superimposition method instead of just (l, m) = (2, 2) dominant mode
# Render needs to resolve sharpness at certain angles and failure to rotate wave propogation

# change directories to whatever
input_directory = r'C:\Users\sethw\OneDrive\PythonFiles\BHVisResearch\scripts\data\S-imposed-strain'
output_directory = "./new/outputmeshstates/"

# Reads inputted strain data - returns data file row length, initial strain data (complex), and time values

def initialize(input_path):
    with open(input_path, "r") as f:
        strain_data = np.array([list(map(float, line.split())) for line in f if not(line.startswith("#"))])

    strain_data = np.unique(strain_data, axis=0)
    length = len(strain_data)
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j * h_imag
    if not (np.diff(h_time) > 0).all():
        raise ValueError("h_time must be strictly increasing")
    return length, np.array(h_strain), np.array(h_time)

def main():
    # File parameters
    input_file = "./r100/Rpsi4_l2-r0100.0.txt"
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

    # loop/read through strain files and store
    #input_file_pattern = r"S-imposed_strain_theta(\d+)-r0100.0.txt"
    input_file_pattern = re.compile(r'theta(\d+)')
    data_dict = {}
    for current_file in os.listdir(input_directory):    
        #match = re.search(input_file_pattern, current_file)
        match = input_file_pattern.search(current_file)
        if match:
            theta = int(match.group(1))
            file_path = os.path.join(input_directory, current_file)
            data_dict[theta] = initialize(file_path)
            # access values through data_dict[theta][0 for length, 1 for strain, 2 for times]
    print("storage done")
    
    # Plot parameters
    numRadius = 450  # Number of points in the radial direction
    states_per_angle = data_dict[0][0] # number of time states in each saved strain angle
    print(type(states_per_angle))
    numTheta = len(data_dict)  # Number of points per circle
    display_radius = 300  # Mesh radius
    R_ext = 100  # Extraction radius
    s = -2
    l = 2
    m = 2
    scale_factor = 600
    omit_radius = 4
    ell_max = 8
    
    percentage = np.round(np.linspace(0, states_per_angle, 101)).astype(int)

    # Pre-compute the theta and radius values
    radius_values = np.linspace(0, display_radius, numRadius)
    azimuth_values = np.linspace(0, 2 * np.pi, numTheta, endpoint=False)
    colatitude = np.pi/2
    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)
    print("interpolation begin")
    #generate and filter target times to pre-interpolate
    #h_array = np.zeros(states_per_angle) # stores time interpolated strain for each theta point as index
    
    length = data_dict[0][0]
    h_time = data_dict[0][2]
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
    
    h_array = {}
    for theta_pt, values_tuple in data_dict.items():
        length, h_strain, h_time, = values_tuple
        # h_array[theta_point][time]
        h_array[theta_pt] = np.interp(t_array, h_time, h_strain)
    '''
    NOTE: h_array now stores the strains in a dictionary, structured [theta point][time state][radius]
            in the future we'll need a fourth dimension for phi (polar values)
    '''
    print("interpolation done")
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
            for i, azi in enumerate(azimuth_values): # i corresponds to theta_pt
                h_tR = h_array[i][state][j]
                x = x_values[j, i]
                y = y_values[j, i]
                if radius <= omit_radius: # Introduce a discontinuity to make room for the Black Holes
                    z = np.nan
                    strain_value = np.nan
                else:
                    strain_value = h_tR.real
                    z = strain_value * scale_factor
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