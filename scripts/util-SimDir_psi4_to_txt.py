# Created to plot specific modes of Psi4 and save data to txt file
#
# Note: This script is designed to work with ET simulation directories

import sys
import matplotlib.pyplot as plt
import numpy as np
from kuibit.simdir import SimDir as sd

# Check if the right number of command-line arguments is provided
if len(sys.argv) < 5:
    print('Usage: python script.py <simulation_directory> <l_value> <m_value> <save_file>')
    sys.exit(1)

# The first command-line argument is the simulation directory
sim_dir = sys.argv[1]

# The second and third command-line arguments are l and m
l = int(sys.argv[2])
m = int(sys.argv[3])

# The fourth command-line argument is whether to save the file
save_file = sys.argv[4].lower() == 'true'

sim = sd(sim_dir)

gws = sim.gws

r = 100

# Use l and m here
psi4_lm = gws[r][(l, m)]

# If save_file is True, save the file
if save_file:
    filename = f"psi4_{l}{m}.txt"
    with open(filename, "w") as f:
        f.write("# Time_Domain_Psi4 Real Imaginary\n")
        np.savetxt(f, np.column_stack((psi4_lm.t, psi4_lm.values.real, psi4_lm.values.imag)), fmt='%-20.15f', delimiter=' ')

# New plotting section
plt.figure(figsize=(10, 6))
plt.plot(psi4_lm.t, psi4_lm.values.real)
plt.title(f'Real Part of Psi4 ({l},{m}) Mode')
plt.xlabel('Time')
plt.ylabel('Real(Psi4)')
plt.grid(True)
plt.show()