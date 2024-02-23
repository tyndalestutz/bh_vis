import scipy 
from math import pi,cos,sin,log
import matplotlib.pyplot as plt


#error function
def ERF(x,x0,w):
    return 0.5 * (scipy.special.erf(((x)-(x0))/(w)) + 1.0)

t_0 = 0
t_final = 2000
num_data_points = 6578
deltat = (t_final - t_0)/num_data_points
initial_seperation = 10
orbital_period = 225
omega = 2*pi/orbital_period

bh_data = [[],[],[],[]]

with open("bh_synthetic.csv", "w") as file:
    outstr = "time,BH1x,BH1y,BH1z,BH2x,BH2y,BH2z\n"
    file.write(outstr)
    for i in range(num_data_points):
        time = t_0 + deltat * i 
        orbital_separation = ERF(time, 1000, -10)
        #BH1 coords
        BH1x = initial_seperation/2 * orbital_separation * cos(omega * time) * (1 - i/num_data_points)
        BH1y = initial_seperation/2 * orbital_separation * sin(omega * time) * (1 - i/num_data_points) 
        BH1z = 0
        #BH2 coords
        BH2x = -initial_seperation/2 * orbital_separation * cos(omega * time) * (1 - i/num_data_points)
        BH2y = -initial_seperation/2 * orbital_separation * sin(omega * time) * (1 - i/num_data_points)
        BH2z = 0
        bh_data[0].append(BH1x)
        bh_data[1].append(BH1y)
        bh_data[2].append(BH2x)
        bh_data[3].append(BH2y)
        #typecast
        outstr = str(time) + "," + str(BH1x) + "," + str(BH1y) + "," + str(BH1z) + "," + str(BH2x) + "," + str(BH2y) + "," + str(BH2z) + "\n"
        file.write(outstr)

plt.figure(figsize=(10, 6))
plt.scatter(bh_data[0], bh_data[1], label="Black Hole 1")
plt.scatter(bh_data[2], bh_data[3], label="Black Hole 2")
plt.xlabel("X-position")
plt.ylabel("Y-position")
plt.title("Black Hole Trajectories")
plt.legend()
plt.grid(True)

plt.show()
