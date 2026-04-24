'''Summarize all SUMMA Actors logs in a folder.
Summary file is placed inside the log folder. Specifying a summary file name is optional.
Usage: python actors_summarize_logs.py [log_folder] [name_of_summary_file.txt]'''

import os
import re
import sys
import statistics as sts

# ----------------------
# Set defaults
summaryFile = '_actors_log_summary.txt' # default, placed at the top of the log folder

# Handle input arguments
if len(sys.argv) == 1: # sys.argv only contains the script name
    sys.exit('Error: no input folder specified')
else:
    log_folder = sys.argv[1]
    if len(sys.argv) >= 3:
        summaryFile = sys.argv[2]

output_file = os.path.join(log_folder, summaryFile)
# End of input arguments
# ----------------------

# Define the patterns to search for
duration_pattern = re.compile(r'Total Duration = ([\d.]+) Hours')
failed_pattern = re.compile(r'Num Failed = (\d+)')
restarts_pattern = re.compile(r'Num Restarts = (\d+)')
file_pattern = re.compile(r'File Manager Path: .*/fileManager_([^/]+)\.txt')
gru_pattern = re.compile(r'Starting SUMMA Actor, start_gru (\d+), num_gru (\d+)')

# Initialize lists to store the extracted values
durations = []
failures = []
restarts = []
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
            restarts_match = restarts_pattern.search(content)
            file_match = file_pattern.search(content)
            gru_match = gru_pattern.search(content)
            if duration_match and failed_match and restarts_match and file_match and gru_match:
                durations.append(duration_match.group(1))
                failures.append(failed_match.group(1))
                restarts.append(restarts_match.group(1))
                file_names.append(file_match.group(1))
                start_grus.append(gru_match.group(1))
                num_grus.append(gru_match.group(2))
                start_gru = int(gru_match.group(1))
                num_gru = int(gru_match.group(2))
                array = (start_gru - 1) / num_gru
                results.append((filename, f'{duration_match.group(1):>10} Hours', f'{failed_match.group(1):>3} failures', f'{restarts_match.group(1):>3} restarts', f'file {file_match.group(1):<8}', f'array {array:.2f}', f'start_gru {start_gru:>6}', f'num_gru {num_gru:>5}'))
            elif file_match and gru_match:
                # Print the filename if the patterns were not found, didn't finish
                start_gru = int(gru_match.group(1))
                num_gru = int(gru_match.group(2))
                array = (start_gru - 1) / num_gru
                results.append((filename, 'NOT FINISHED', '', f'file {file_match.group(1):<8}', f'array {array:.2f}', f'start_gru {start_gru:>6}', f'num_gru {num_gru:>5}'))
            else:
                # Print the filename if the patterns were not found
                results.append((filename, 'NOT FOUND', '', '', '', '', ''))

# Sort the results by file_match, then array before printing
results.sort(key=lambda x: (x[3], x[4]))

# Calculate summary statistics
total_success = sum(1 for result in results if 'Hours' in result[1])
total_summa = sum(1 for result in results if result[1] == 'NOT FINISHED')
total_other = sum(1 for result in results if result[1] == 'NOT FOUND')
total = total_success + total_summa + total_other

if total > 0:
    pct_success = total_success / total * 100
    pct_summa = total_summa / total * 100
    pct_other = total_other / total * 100
else:
    pct_success = 0
    pct_summa = 0
    pct_other = 0

# Calculate computation time stats
computation_time = [float(time) for time in durations if float(time) >= 0]

if computation_time:
    st_min = min(computation_time)
    st_max = max(computation_time)
    st_mean = sts.mean(computation_time)
    st_median = sts.median(computation_time)
else:
    st_min = None
    st_max = None
    st_mean = None
    st_median = None

# Print the sorted results
for result in results:
    print(f'{result[0]:<15} | {result[1]:>10} | {result[2]:>3} | {result[3]:>3} | {result[4]:<8} | {result[5]:<8} | {result[6]:<12} | {result[7]:<8}')

    print('\nSuccess stats')
    print('Success' + '\t \t \t \t {:.2f}% \n'.format(pct_success))
    print('SUMMA error' + '\t \t \t {:.2f}% \n'.format(pct_summa))
    print('Early termination' + '\t {:.2f}% \n'.format(pct_other))
    print('\nTime needed for successful computations')
    if computation_time:        
        print('Min time ' + '\t \t \t {:.2f} h \n'.format(st_min))
        print('Median time ' + '\t \t {:.2f} h \n'.format(st_median))
        print('Mean time ' + '\t \t \t {:.2f} h \n'.format(st_mean))
        print('Max time ' + '\t \t \t {:.2f} h \n'.format(st_max))
    else:
        print('No successful computations found\n')

# Write the sorted results to the output file
with open(output_file, 'w') as file:
    for result in results:
        file.write(f'{result[0]:<15} | {result[1]:>10} | {result[2]:>3} | {result[3]:>3} | {result[4]:<8} | {result[5]:<8} | {result[6]:<12} | {result[7]:<8}\n')

    file.write('\nSuccess stats\n')
    file.write('Success' + '\t \t \t \t {:.2f}% \n'.format(pct_success))
    file.write('SUMMA error' + '\t \t \t {:.2f}% \n'.format(pct_summa))
    file.write('Early termination' + '\t {:.2f}% \n'.format(pct_other))
    file.write('\nTime needed for successful computations\n')
    if computation_time:
        file.write('Min time ' + '\t \t \t {:.2f} h \n'.format(st_min))
        file.write('Median time ' + '\t \t {:.2f} h \n'.format(st_median))
        file.write('Mean time ' + '\t \t \t {:.2f} h \n'.format(st_mean))
        file.write('Max time ' + '\t \t \t {:.2f} h \n'.format(st_max))
    else:
        file.write('No successful computations found\n')

print(f'Summary written to {output_file}')