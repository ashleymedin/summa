# regroup summa runs output into different group size
# written by A. Van Beusekom (2025)

# Run:
# python regroup_summa.py sundials_1en8

import os
from glob import glob
import netCDF4 as nc
import numpy as np
import sys
import multiprocessing as mp

# number of groups to divide the output into
num_groups = 65
chunk_size = 1000  # Number of records to process at a time
gru_num = 517315  # Number of GRUs 

missing = False  # if appending nan hrus to batch because failed
missgru = 72055933  # batch 205 summa-be32 value
misshru = missgru  # could be different

run_local = False
if run_local:
    top_fold = '/Users/amedin/Research/USask/test_py/'
    method_name = 'sundials_1en8'
else:
    top_fold = '/project/k/kshook/avanb/enthalpy_paper/runs/'
    method_name = sys.argv[1]  # sys.argv values are strings by default so this is fine (sundials_1en8 or be64)

ncdir = top_fold + 'summa-' + method_name + '_nocat'
file_pattern = 'run1_G*_timestep.nc'
ctdir = top_fold + 'summa-' + method_name

# Get list of split summa output files (hardwired pattern)
outfilelist0 = glob((ncdir + '/' + file_pattern))
outfilelist0.sort()

# -- functions
def concatenate_files_in_range(outfilelist0, ctdir, start_gru, end_gru):
    out_name = f'all_G{start_gru + 1:06d}-{end_gru:06d}.nc'
    gru_num = 0
    hru_num = 0

    # Filter files within the specified range
    filtered_files = [file for file in outfilelist0 if start_gru <= int(file.split('_')[1][1:7]) <= end_gru]

    for file in filtered_files:
        f = nc.Dataset(file)
        gru_num += len(f.dimensions['gru'])
        hru_num += len(f.dimensions['hru'])

    # Write output
    with nc.Dataset(filtered_files[0]) as src:
        with nc.Dataset(ctdir + '/' + out_name, "w") as dst:
            # Copy dimensions
            for name, dimension in src.dimensions.items():
                if name == 'gru':
                    dst.createDimension(name, gru_num)
                elif name == 'hru':
                    dst.createDimension(name, hru_num)
                else:
                    dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            # Copy variable attributes all at once via dictionary
            gru_vars = []  # variable name, gru axis in variable dimension for concatenation.
            hru_vars = []
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                # Note here the variable dimension name is the same, but size has been updated for gru and hru.

                # Assign different values depending on dimension
                dims = variable.dimensions
                if 'gru' in dims:
                    gru_vars.append([name, dims.index('gru')])
                elif 'hru' in dims:
                    hru_vars.append([name, dims.index('hru')])
                else:
                    dst[name][:] = src[name][:]

            gru_vars_num = len(gru_vars)
            hru_vars_num = len(hru_vars)

            # Read and write values for gru and hru dimensioned variables in chunks
            for i, file in enumerate(filtered_files):
                print("combining file %d %s" % (i, file))
                f = nc.Dataset(file)
                for j in range(gru_vars_num):
                    gru_var_name = gru_vars[j][0]
                    for start in range(0, len(f[gru_var_name]), chunk_size):
                        end = min(start + chunk_size, len(f[gru_var_name]))
                        dst.variables[gru_var_name][start:end] = f[gru_var_name][start:end]

                for j in range(hru_vars_num):
                    hru_var_name = hru_vars[j][0]
                    for start in range(0, len(f[hru_var_name]), chunk_size):
                        end = min(start + chunk_size, len(f[hru_var_name]))
                        dst.variables[hru_var_name][start:end] = f[hru_var_name][start:end]

            # If missing HRUs, this is slow or broken
            if missing:
                new_index = np.append(dst["gru"][:], missgru)
                dst.reindex({"gru": new_index})
                dst.sel(gru=missgru)["gruId"] = missgru

                new_index = np.append(dst["hru"][:], misshru)
                dst.reindex({"hru": new_index})
                dst.sel(gru=misshru)["hruId"] = misshru

        print("wrote output: %s" % (ctdir + '/' + out_name))

    return out_name

def divide_into_groups(all_file, ctdir, start_gru, end_gru):
    with nc.Dataset(ctdir + '/' + all_file) as src:
        hru_num = gru_num
        start_hru = start_gru
        end_hru = end_gru
        out_name = f'run1__G{start_gru + 1:06d}-{end_gru:06d}_timestep.nc'

        with nc.Dataset(ctdir + '/' + out_name, "w") as dst:
            # Copy dimensions
            for name, dimension in src.dimensions.items():
                if name == 'gru':
                    dst.createDimension(name, end_gru - start_gru)
                elif name == 'hru':
                    dst.createDimension(name, end_hru - start_hru)
                else:
                    dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            # Copy variable attributes all at once via dictionary
            for name, variable in src.variables.items():
                dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)

                # Assign different values depending on dimension
                dims = variable.dimensions
                if 'gru' in dims:
                    for start in range(start_gru, end_gru, chunk_size):
                        end = min(start + chunk_size, end_gru)
                        dst[name][start - start_gru:end - start_gru] = variable[start:end]
                elif 'hru' in dims:
                    for start in range(start_hru, end_hru, chunk_size):
                        end = min(start + chunk_size, end_hru)
                        dst[name][start - start_hru:end - start_hru] = variable[start:end]
                else:
                    dst[name][:] = variable[:]

            print(f"wrote output: {ctdir + '/' + out_name}")
            # remove the all file
            os.remove(ctdir + '/' + all_file)

# -- end functions

def process_group(g):
    size_g = int(np.ceil(gru_num / num_groups))
    start_gru = g * size_g
    end_gru = min((g + 1) * size_g, gru_num)
    all_file = concatenate_files_in_range(outfilelist0, ctdir, start_gru, end_gru)
    divide_into_groups(all_file, ctdir, start_gru, end_gru)

# -- end functions

if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(process_group, range(num_groups))