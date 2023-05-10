import vtk, os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm

'''
To Do:
- Render whole strain file
- Compile all l, m, strain values given
- Add functionality to customize output directory
'''


# File parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/bh_vis/scripts/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
#output_directory = parent_directory + "/../../data/mesh/" + folder_name
output_directory = "/home/guest/Documents/BH_Vis_local/data/mesh/gw_test/" #i know this is cheating


# Set up grid parameters
status_messages = True
plot_strain = True
l = 2
m = 2
num_points_x = 200
num_points_y = 200
num_points_z = 1
R_ext = 100 # Radius of data extraction we are taking the gw strain from

# Calculate spin-weighted spherical harmonics for every point in the mesh, returns them as complex values in a 2d array
def set_sph_harm_array(l,m,s):
    global num_points_x, num_points_y

    # Store values in 2d array
    sph_harm_points = np.zeros((num_points_x,num_points_y),dtype=np.complex128) 

    for j in range(num_points_y):
        y = -(num_points_y - 1) + 2 * j
        for i in range(num_points_x):
            x = -(num_points_x - 1) + 2 * i
            r = np.sqrt(x**2 + y**2)
            theta = np.pi/2 # For 2d purposes, theta (polar angle) will be constant
            phi = np.arctan2(y, x)
            Y_lm = sph_harm(l, m, phi, theta) # Calculate spherical harmonics using scipy.special funciton
            Y = (-1)**s * np.sqrt((2*l+1)/(4*np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m)) * Y_lm # incorporate spin-weight factor
            sph_harm_points[i,j] = Y
    return sph_harm_points


# This function allows a linearly interpolated strain to be found given a target time
def interpolated_strain(target_time, source_time, data):
    # Check for data errors before interpolating
    if len(source_time) != len(data):
        raise ValueError("source_time and data must have the same number of rows")
    if not (np.diff(source_time) > 0).all():
        raise ValueError("source_time must be strictly increasing")
    
    # Interpolate the data using numpy's interp function
    interpolated_data = np.interp(target_time, source_time, data)
    return interpolated_data

# Reads inputted strain data - returns data file row length, initial strain data (complex), spin-weighted spherical harmonics (complex), and time values
def initialize():
    # If output directory exists and has items, asks before the program overwrites. If the directory does not exist, the program exits
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
        
    # Function to ignore lines starting with #
    def valid_line(line):
        return not line.startswith("#")
    
    # Load the psi4 strain data into arrays - already precalculated
    with open(input_file, 'r') as f:
        strain_data = np.array([list(map(float, line.split())) for line in f if valid_line(line)])
    
    # Sort by time and remove duplicates
    strain_data = np.unique(strain_data, axis=0)

    # Separate time, real and imainary parts of h
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j*h_imag
    

    length = len(h_strain)

    # Set spin-weighted spherical harmonics
    sph_harm_points = set_sph_harm_array(2,2,-2)

    return length, h_strain, sph_harm_points, h_time

# Assign values from initialize() function
length, h_strain, sph_harm_points, h_time = initialize()

# Plot the strain in 2d without spin-weighed spherical harmonics if wanted
def show_strain_plot():
    plt.plot(h_time, h_strain.real)
    plt.title("Gravitational Wave Strain vs Time\nl =" + str(l) + "m =" + str(m) + "\nData Extraction Radius =" + str(R_ext) + "meters")
    plt.xlabel("Time")
    plt.ylabel("Real Part of Strain")

    plt.show()
if plot_strain:
    show_strain_plot()

start_time = time.time() # Start timer
percentage = np.round(np.linspace(0, length, 101)).astype(int) #creates an array of 10% points

# Loop through all time values in the data
state = 0 # For file naming purposes (***probably a better way to name them than this)
for current_time in h_time:
    state += 1
    # Create mesh
    points = vtk.vtkPoints()
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(num_points_x, num_points_y, num_points_z)
    
    # Output data generation progress to terminal - currently broken
    t = np.where(h_time == current_time)[0][0]
    if status_messages and t == 10:
        end_time = time.time() #end timer 
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}")
    if status_messages and t != 0 and np.isin(t,percentage):
        print(f" {int(t * 100 / (length - 1))}% done", end="\r") #create percentage status message ### THIS ISNT WORKING????
    
    # For every time value, set up x, y mesh points
    for j in range(num_points_y):
        y = -(num_points_y - 1) + 2 * j
        for i in range(num_points_x):
            x = -(num_points_x - 1) + 2 * i
            current_r = np.sqrt(x**2 + y**2)
            # For every point in the mesh, calculate the adjusted time to find the strain at based on radius, simulation time, and extraction radius
            target_time = current_time - current_r + R_ext
            
            # Find initial and final times in the data, constrain target_time within those values
            time_0 = np.min(h_time)
            time_f = np.max(h_time)
            if (target_time < time_0):
                target_time = time_0
            elif (target_time > time_f):
                target_time = time_f

            # Find the intermediate strain that's dependent on current time and data extraction radius
            h_tR = interpolated_strain(target_time, h_time, h_strain)
            # Find the spin weighted spherical harmonic from the data array made for the current point
            Y = sph_harm_points[i, j]

            #Plot z based on the real part of the product of h_tR and Y
            z = Y.real*h_tR.real - Y.imag*h_tR.imag
            points.InsertNextPoint(x, y, z*100)
            
    grid.SetPoints(points)

    # Write mesh to file
    writer = vtk.vtkXMLStructuredGridWriter()
    filename = output_directory + f"/state{state}.vts"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()

print("Mesh database completed in",output_directory)