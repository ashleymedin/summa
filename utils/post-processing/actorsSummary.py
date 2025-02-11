import os
import re

# Define the folder containing the log files
log_folder = os.path.expanduser('~/summaNorthAmerica_settings/logs')
output_file = 'summary.txt'

# Define the patterns to search for
duration_pattern = re.compile(r'Total Duration = ([\d.]+) Hours')
failed_pattern = re.compile(r'Num Failed = (\d+)')
file_pattern = re.compile(r'File Manager Path: /home/x-avanb/summaNorthAmerica_settings/fileManager_([^/]+)\.txt')
gru_pattern = re.compile(r'Starting SUMMA Actor, start_gru (\d+), num_gru (\d+)')

# Initialize lists to store the extracted values
durations = []
failures = []
file_names = []
start_grus = []
num_grus = []

# Iterate over all files in the folder
for filename in os.listdir(log_folder):
    if filename.startswith('log'):
        filepath = os.path.join(log_folder, filename)
        with open(filepath, 'r') as file:
            content = file.read()
            # Search for the patterns in the file content
            duration_match = duration_pattern.search(content)
            failed_match = failed_pattern.search(content)
            file_match = file_pattern.search(content)
            gru_match = gru_pattern.search(content)
            if duration_match and failed_match and file_match and gru_match:
                durations.append(duration_match.group(1))
                failures.append(failed_match.group(1))
                file_names.append(file_match.group(1))
                start_grus.append(gru_match.group(1))
                num_grus.append(gru_match.group(2))
                start_gru = int(gru_match.group(1))
                num_gru = int(gru_match.group(2))
                array = (start_gru - 1) / num_gru
                print(f'{filename:<15} | {duration_match.group(1):>10} Hours | {failed_match.group(1):>3} failures | file {file_match.group(1):<8} | array {array:.2f} | start_gru {start_gru:>6} | num_gru {num_gru:>5}')
            elif file_match and gru_match:
                # Print the filename if the patterns were not found, didn't finish
                start_gru = int(gru_match.group(1))
                num_gru = int(gru_match.group(2))
                array = (start_gru - 1) / num_gru
                print(f'{filename:<15} | {"NOT FINISHED":>31} | file {file_match.group(1):<8} | array {array:.2f} | start_gru {start_gru:>6} | num_gru {num_gru:>5}')
            else:
                # Print the filename if the patterns were not found
                print(f'{filename:<15} | {"NOT FOUND":>31} |')

# Write the extracted values to the output file
#with open(output_file, 'w') as file:
#    for duration, failure, file_name, start_gru, num_gru in zip(durations, failures, file_names, start_grus, num_grus):
#        file.write(f'Total Duration = {duration} Hours\n')
#        file.write(f'Num Failed = {failure}\n')
#        file.write(f'Created output file: /anvil/scratch/x-avanb/{file_name}\n')
#        file.write(f'Starting SUMMA Actor, start_gru {start_gru}, num_gru {num_gru}\n')
#
#print(f'Summary written to {output_file}')