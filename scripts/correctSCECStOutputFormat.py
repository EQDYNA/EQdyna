import os
import re

# Function to detect and correct invalid scientific notation
def correct_scientific_notation_line(line):
    # Regular expression to detect numbers like 0.1341312-114 or -0.7359989-113
    corrected_line = re.sub(r'(\d*\.\d+)([-+]\d+)', r'\1E\2', line)
    return corrected_line

# Function to process and overwrite the file with corrected content
def process_file(file_path):
    # Temporary output file to store corrected content
    temp_output_file = file_path + '.tmp'
    
    with open(file_path, 'r') as infile, open(temp_output_file, 'w') as outfile:
        for line in infile:
            # Only process lines with data (ignore comment lines)
            if not line.startswith('#'):
                corrected_line = correct_scientific_notation_line(line)
                outfile.write(corrected_line)
            else:
                # Write comment lines as-is
                outfile.write(line)
    
    # Replace the original file with the corrected file
    os.replace(temp_output_file, file_path)
    print(f"Processed and updated file: {file_path}")

# Function to find and process all files starting with 'body' or 'faultst'
def process_all_files_in_directory(directory):
    # Get all files in the directory
    for file_name in os.listdir(directory):
        # Check if the file name starts with 'body' or 'faultst'
        if file_name.startswith('body') or file_name.startswith('faultst'):
            file_path = os.path.join(directory, file_name)
            process_file(file_path)

# Directory where the files are located
directory_path = '.'  # Replace with your actual directory path

# Process all matching files in the directory
process_all_files_in_directory(directory_path)
