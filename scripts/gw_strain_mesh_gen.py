import vtk, os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm


#file parameters
folder_name = "gw_test"
input_file = "/home/guest/Documents/BH_Vis/data/gw_data/Rpsi4_l2-r0100.0_strain.txt"
parent_directory = os.path.dirname(os.path.dirname(__file__))
#output_directory = parent_directory + "/../../data/mesh/" + folder_name
output_directory = "/home/guest/Documents/BH_Vis/data/mesh/gw_test/" #i know this is cheating


# Set up grid parameters
status_messages = True 
l = 2
m = 2
num_points_x = 200
num_points_y = 200
num_points_z = 1
R_ext = 100 # radius of data extraction we are taking the gw strain from

# Most likely this will not be used.
propogation_speed = 1 #how many pixels the gw travels in one delta t, i.e. (c/ units of time per second) * (pixels/meter)

#note: maybe make amplitude part of the render script instead of this script?
amplitude = 100  # Set the amplitude of the wave visualization.

#interpolated array allows for non-integer indecies and negative indecies return 0
######## this is currently unused for most recent render.
class InterpolatedArray:
    def __init__(self, data):
        self.data = data
        
    def __getitem__(self, key):
        data = self.data
        if key > len(data) - 1:
            print(f"IndexError: index {key} is out of bounds for array with size {len(data)}")
            exit()
        elif key < 0: 
            return 0
        elif key == int(key):
            return data[key]
        else:
            slope = (data[int(key) + 1] - data[int(key)])
            return (slope * (key - int(key))) + data[int(key)] #linear interpolation m(x1 - x0) + x0

##WARNING: figure out spin-weoight factor
def set_sph_harm_array(l,m,s):
    global num_points_x, num_points_y

    sph_harm_points = np.zeros((num_points_x,num_points_y),dtype=np.complex128) 

    for j in range(num_points_y):
        y = -(num_points_y - 1) + 2 * j
        for i in range(num_points_x):
            x = -(num_points_x - 1) + 2 * i
            r = np.sqrt(x**2 + y**2)
            theta = np.pi/2
            #theta = np.arccos(2/r) # Etienne said to use this for theta?
            phi = np.arctan2(y, x)
            Y_lm = sph_harm(l, m, phi, theta)
            # compute the spin-weighted spherical harmonic
            Y = (-1)**s * np.sqrt((2*l+1)/(4*np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m)) * Y_lm

            sph_harm_points[i,j] = Y
    print(sph_harm_points)
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

    # separate time, real and imainary parts of h
    h_time, h_real, h_imag = strain_data[:, 0], strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j*h_imag
    

    length = len(h_strain)
    #h_strain = InterpolatedArray(h_strain) #allows non-integer indecies of h_strain
    sph_harm_points = set_sph_harm_array(2,2,-2)

    return length, h_strain, sph_harm_points, h_time
    
length, h_strain, sph_harm_points, h_time = initialize()

start_time = time.time() #start timer
percentage = np.round(np.linspace(0, length, 101)).astype(int) #creates an array of 10% points

def total_strain_real(strain, R_ext, t, x, y):
    r = np.sqrt(x**2 + y**2)
    strain_real = sph_harm_points.real*strain[t].real - sph_harm_points.imag*strain[t].imag
    return strain_real

# loop through all time values in the data
state = 0 # just for file naming purposes
for current_time in h_time:
    state += 1
    # Create mesh
    points = vtk.vtkPoints()
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(num_points_x, num_points_y, num_points_z)
    
    '''
    if status_messages and t == 10:
        end_time = time.time() #end timer 
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}")
    if status_messages and t != 0 and np.isin(t,percentage):
        print(f" {int(t * 100 / (length - 1))}% done", end="\r") #create percentage status message ### THIS ISNT WORKING????
    '''
    # for every time value, set up x, y mesh points
    for j in range(num_points_y):
        y = -(num_points_y - 1) + 2 * j
        for i in range(num_points_x):
            x = -(num_points_x - 1) + 2 * i
            current_r = np.sqrt(x**2 + y**2)
            # for every point in the mesh, calculate the adjusted time to find the strain at based on radius, simulation time, and extraction radius
            target_time = current_time - current_r + R_ext
            
            # find initial and final times in the data, constrain target_time within those values
            time_0 = np.min(h_time)
            time_f = np.max(h_time)
            if (target_time < time_0):
                target_time = time_0
            elif (target_time > time_f):
                target_time = time_f

            # find the intermediate strain that's dependent on current time and data extraction radius
            h_tR = interpolated_strain(target_time, h_time, h_strain)
            # find the spin weighted spherical harmonic from the data array made for the current point
            Y = sph_harm_points[i, j]

            #Plot z based on the real part of the product of h_tR and Y
            z = Y.real*h_tR.real - Y.imag*h_tR.imag
            points.InsertNextPoint(x, y, z*100)
            #print(x, y, z)
            '''
            for time in h_time:

                z = total_strain_real(h_strain, R_ext, time, x, y) * amplitude
                #h_total = sph_harm_points[i,j]*h_strain[t - int(r/propogation_speed)]

                # plot the real part of the strain into z
                #z = h_total.real * amplitude 
                points.InsertNextPoint(x, y, z)
            '''
            
    grid.SetPoints(points)

    # Write mesh to file
    writer = vtk.vtkXMLStructuredGridWriter()
    filename = output_directory + f"/state{state}.vts"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()
print("Mesh database completed in",output_directory)