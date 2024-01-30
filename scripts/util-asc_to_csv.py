# Created to extract specific columns of asc file into csv.
# This script is specifically helpful for converting puncture tracker
# data into csv format, which can then be used with VisIt.

import sys
import csv

def parse_file(input_file, output_file, column_indexes_to_keep):
    headers = []
    data = []

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith("# data columns:"):
                all_headers = line.strip().split()[3:]
                # Only keep the headers for the columns we're interested in
                headers = [all_headers[i] for i in column_indexes_to_keep if i < len(all_headers)]
            elif not line.startswith("#"):
                all_data = line.strip().split()
                # Only keep the data for the columns we're interested in
                row = [all_data[i] for i in column_indexes_to_keep if i < len(all_data)]
                data.append(row)

    with open(output_file, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        writer.writerows(data)


# Check if the right number of command-line arguments is provided
if len(sys.argv) < 3:
    print('Usage: python script.py <input_file> <column_indexes_to_keep>')
    print('Example: python3 asc_to_csv.py combined.asc "8,22,23,32,33"')
    sys.exit(1)

# The first command-line argument is the input filename
input_file = sys.argv[1]

# The second command-line argument is the string of column indexes to keep
column_indexes_to_keep = list(map(int, sys.argv[2].split(',')))

# Generate the output filename based on the input filename
# If the input filename is 'combined.asc',
# the output filename will be 'combined_converted1.csv'
output_file = input_file.rsplit('.', 1)[0] + '.csv'

parse_file(input_file, output_file, column_indexes_to_keep)