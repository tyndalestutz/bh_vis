# Initialize an empty dictionary to store the strain data
strain_data = {}

# Loop over all values of l 
for l in range(2, 9):
    # Construct the filename based on l
    filename = f"path_to_strain_files/Rpsi4_l{l}_strain.txt"
    
    # Load the data from the file
    with open(filename, 'r') as f:
        data = np.array([list(map(float, line.split())) for line in f if not line.startswith("#")])

    # Extract the time column
    time = data[:, 0]
    
    # Loop over all values of m and extract corresponding strain data
    for m in range(-l, l+1):
        # Calculate the index of the mth mode in the data array
        m_index = 2 * (m + l) + 1  # Add 1 because the first column is time
        
        # Extract the real and imaginary parts for the mth mode
        real_part = data[:, m_index]
        imag_part = data[:, m_index + 1]

        # Combine the real and imaginary parts to get the complex strain
        strain = real_part + 1j * imag_part

        # Store the strain data in the dictionary
        strain_data[(l, m)] = (time, strain)
