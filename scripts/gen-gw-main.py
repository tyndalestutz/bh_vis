import os, time, vtk, numpy as np, quaternionic, spherical

# Reads inputted strain data - returns data file row length, initial strain data (complex), and time values
def initialize(input_file, output_directory):
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

def main():
    # File parameters
    input_file = "./Rpsi4_l2-r0100.0_strain.txt"
    output_directory = "./new/outputmeshstates-l2m2/"
    status_messages = True

    length, h_strain, h_time = initialize(input_file, output_directory)

    # Plot parameters
    numRadius = 450  # Number of points in the radial direction
    numTheta = 180  # Number of points per circle
    display_radius = 300  # Mesh radius
    R_ext = 100  # Extraction radius
    s = -2
    l = 2
    m = 2
    scale_factor = 600
    omit_radius = 4
    ell_max = 8
    D = spherical.Wigner(ell_max)
    
    percentage = np.round(np.linspace(0, length, 101)).astype(int)

    # Pre-compute the theta and radius values
    radius_values = np.linspace(0, display_radius, numRadius)
    azimuth_values = np.linspace(0, 2 * np.pi, numTheta, endpoint=False)
    colatitude = np.pi/2
    rv, az = np.meshgrid(radius_values, azimuth_values, indexing="ij")
    x_values = rv * np.cos(az)
    y_values = rv * np.sin(az)
    sY =  D.sYlm(s, quaternionic.array.from_spherical_coordinates(colatitude, azimuth_values))

    #generate and filter target times to pre-interpolate #Revisit
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
            for i, azi in enumerate(azimuth_values):
                x = x_values[j, i]
                y = y_values[j, i]
                swsh = sY[i][D.Yindex(l,m)]
                if radius <= omit_radius: # Introduce a discontinuity to make room for the Black Holes
                    z = np.nan
                    strain_value = np.nan
                else:
                    strain_value = swsh.real * h_tR.real - swsh.imag * h_tR.imag #ignores imaginary cross-terms
                    #strain_value = (swsh * h_tR).real
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