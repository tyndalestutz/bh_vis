from scipy import special
from math import pi,cos,sin
import os

parent_directory = os.path.dirname(os.path.dirname(__file__))

#error function
def ERF(x,x0,w):
    return 0.5 * (special.erf(((x)-(x0))/(w)) + 1.0)

destination_directory = parent_directory + r"/data/synthetic_coords/"
filename = "synthetic_data_ang_momentum.csv"
t_final = 2000
num_data_pts = 1000
orbital_period = 225
radius = 5
omega = 2*pi/orbital_period

deltat = (t_final)/num_data_pts
with open(destination_directory + filename, 'w') as file:
    file.write("time,BH1x,BH1y,BH1z,L1x,L1y,L1z,BH2x,BH2y,BH2z,L2x,L2y,L2z\n")
    bh1 = [0,0,0] #in the form [x,y,z]
    bh2 = [0,0,0] #in the form [x,y,z]
    L1 = [0,0,0] #in the form [|x|,|y|,|z|]
    L2 = [0,0,0] #in the form [|x|,|y|,|z|]
    for i in range(num_data_pts):
        time = deltat * i 
        orbital_separation = ERF(time, 1000, -500) 
        #BH1 data
        bh1[0] =  radius * orbital_separation * cos(omega * time) 
        bh1[1] = radius * orbital_separation * sin(omega * time)
        bh1[2] = 1
        L1[0] = 0 
        L1[1] = -1
        L1[2] = 1
        #BH2 data
        bh2[0] = -radius * orbital_separation * cos(omega * time)
        bh2[1] = -radius * orbital_separation * sin(omega * time)
        bh2[2] = 1
        L2[0] = 0
        L2[1] = 1.5
        L2[2] = 0

        #typecast
        outstr = str(time) + ","
        for i in [bh1,L1,bh2,L2]:
            for j in i:
                outstr += str(j) + "," 
        outstr += "\n"

        file.write(outstr)