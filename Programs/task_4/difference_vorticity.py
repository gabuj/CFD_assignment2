import numpy as np

def calculate_file_difference(file1_path, file2_path, output_path):
    try:
        # Load the data from both files
        # It handles whitespace automatically and supports scientific notation
        data1 = np.loadtxt(file1_path)
        data2 = np.loadtxt(file2_path)

        # Check if files have the same number of rows
        if data1.shape != data2.shape:
            print(f"Error: The files have different dimensions ({data1.shape} vs {data2.shape}).")
            return

        # Create the result array:
        # Columns 1 and 2 are taken from the first file
        # Column 3 is the difference (File 1 - File 2)
        result = np.zeros_like(data1)
        result[:, 0] = data1[:, 0]
        result[:, 1] = data1[:, 1]
        result[:, 2] = data1[:, 2] - data2[:, 2]

        # Save to a new file
        # The fmt string ensures high precision and scientific notation matching your file
        np.savetxt(output_path, result, fmt='%25.18E %25.18E %25.18E')
        
        print(f"Success! Difference file saved to: {output_path}")

    except Exception as e:
        print(f"An error occurred: {e}")

# Usage
file1 = 'vort_4_4'        # Name of your first file
file2 = 'vort0004'      # Name of your second file
output = 'vort_diff_4'      # Name of the output file

calculate_file_difference(file1, file2, output)