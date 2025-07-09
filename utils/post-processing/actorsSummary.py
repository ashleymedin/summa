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

# Initialize a list to store the results
results = []

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
                results.append((filename, f'{duration_match.group(1):>10} Hours', f'{failed_match.group(1):>3} failures', f'file {file_match.group(1):<8}', f'array {array:.2f}', f'start_gru {start_gru:>6}', f'num_gru {num_gru:>5}'))
            elif file_match and gru_match:
                # Print the filename if the patterns were not found, didn't finish
                start_gru = int(gru_match.group(1))
                num_gru = int(gru_match.group(2))
                array = (start_gru - 1) / num_gru
                results.append((filename, 'NOT FINISHED', '', f'file {file_match.group(1):<8}', f'array {array:.2f}', f'start_gru {start_gru:>6}', f'num_gru {num_gru:>5}'))
            else:
                # Print the filename if the patterns were not found
                results.append((filename, 'NOT FOUND', '', '', '', '', ''))

# Sort the results by logID
results.sort(key=lambda x: x[0])

# Print the sorted results
for result in results:
    print(f'{result[0]:<15} | {result[1]:>10} | {result[2]:>3} | {result[3]:<8} | {result[4]:<8} | {result[5]:<12} | {result[6]:<8}')

# Write the sorted results to the output file
with open(output_file, 'w') as file:
    for result in results:
        file.write(f'{result[0]:<15} | {result[1]:>10} | {result[2]:>3} | {result[3]:<8} | {result[4]:<8} | {result[5]:<12} | {result[6]:<8}\n')

print(f'Summary written to {output_file}')