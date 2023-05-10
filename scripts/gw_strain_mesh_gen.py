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
R_ext = 100 # radius we are taking the gw strain from
propogation_speed = 1 #how many pixels the gw travels in one delta t, i.e. (c/ units of time per second) * (pixels/meter)

#note: maybe make amplitude part of the render script instead of this script?
amplitude = 1000  # Set the amplitude of the wave visualization

#interpolated array allows for non-integer indecies and negative indecies return 0
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
            #theta = np.pi/2 #theta = pi/2
            theta = np.arccos(2/r) # Etienne said to use this for theta
            phi = np.arctan2(y, x)

            if phi < 0:
                phi += 2*np.pi   # Adjust the range to [0, 2pi]
                
            sph_harm_points[i,j] = sph_harm(l,m,theta,phi)
    return sph_harm_points

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
    h_real, h_imag = strain_data[:, 1], strain_data[:, 2]
    h_strain = h_real + 1j*h_imag
    h_time = strain_data[:, 0]


    length = len(h_strain)
    h_strain = InterpolatedArray(h_strain) #allows non-integer indecies of h_strain
    ###### To do: make an interpolate to time function that accepts a time, returns the strain according to that time

    sph_harm_points = set_sph_harm_array(2,2,-2)

    return length, h_strain, sph_harm_points, h_time
    
length, h_strain, sph_harm_points, h_time = initialize()

start_time = time.time() #start timer
percentage = np.round(np.linspace(0, length, 101)).astype(int) #creates an array of 10% points

def total_strain_real(strain, R_ext, t, x, y):
    r = np.sqrt(x**2 + y**2)
    strain_real = sph_harm_points.real*strain[t].real - sph_harm_points.imag*strain[t].imag
    return strain_real

for t, h in enumerate(h_strain):
    # Create mesh
    points = vtk.vtkPoints()
    grid = vtk.vtkStructuredGrid()
    grid.SetDimensions(num_points_x, num_points_y, num_points_z)
    
    if status_messages and t == 10:
        end_time = time.time() #end timer 
        eta = (end_time - start_time) * length / 10
        print(f"Creating {length} meshes and saving them to {output_directory}.\nEstimated time: {eta}")
    if status_messages and t != 0 and np.isin(t,percentage):
        print(f" {int(t * 100 / (length - 1))}% done", end="\r") #create percentage status message ### THIS ISNT WORKING????

    for j in range(num_points_y):
        y = -(num_points_y - 1) + 2 * j
        for i in range(num_points_x):
            x = -(num_points_x - 1) + 2 * i
            r = np.sqrt(x**2 + y**2)

            for time in h_time:

                z = total_strain_real(h_strain, R_ext, time, x, y) * amplitude
                #h_total = sph_harm_points[i,j]*h_strain[t - int(r/propogation_speed)]

                # plot the real part of the strain into z
                #z = h_total.real * amplitude 
                points.InsertNextPoint(x, y, z)
            
    grid.SetPoints(points)

    # Write mesh to file
    writer = vtk.vtkXMLStructuredGridWriter()
    filename = output_directory + f"/state{t-1}.vts"
    writer.SetFileName(filename)
    writer.SetInputData(grid)
    writer.Write()
print("Mesh database completed in",output_directory)