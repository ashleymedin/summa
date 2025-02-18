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
    out_name = f'run1__G{start_gru:06d}-{end_gru:06d}_timestep.nc'
    gru_num = 0
    hru_num = 0

    # Filter files within the specified range
    filtered_files = [file for file in outfilelist0 if start_gru <= int(file.split('_')[4][1:7]) <= end_gru]
    start_gru_old = int(filtered_files[0].split('_')[4][1:7])
    end_gru_old = int(filtered_files[-1].split('_')[4][8:14])

    for file in filtered_files:
        f = nc.Dataset(file)
        gru_num += len(f.dimensions['gru'])
        hru_num += len(f.dimensions['hru'])

    # Write output
    with nc.Dataset(filtered_files[0]) as src:
        with nc.Dataset(ctdir + '/' + out_name, "w") as dst:
            # copy dimensions
            for name, dimension in src.dimensions.items():
                if name == 'gru':
                    dst.createDimension(name, gru_num)
                elif name == 'hru':
                    dst.createDimension(name, hru_num)
                else:
                    dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            # copy variable attributes all at once via dictionary
            gru_vars = [] # variable name, gru axis in variable dimension for concatenation.
            hru_vars = []
            for name, variable in src.variables.items():
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
                dst[name].setncatts(src[name].__dict__)
                # Note here the variable dimension name is the same, but size has been updated for gru and hru.

                # Assign different values depending on dimension
                dims = variable.dimensions
                if 'gru' in dims:
                    gru_vars.append([name,dims.index('gru')])
                elif 'hru' in dims:
                    hru_vars.append([name,dims.index('hru')])
                else:
                    dst[name][:]=src[name][:]

            # read values for gru and hru dimensioned variables
            Dict = {}
            gru_vars_num = len(gru_vars)
            hru_vars_num = len(hru_vars)
            for i, file in enumerate(filtered_files):

                print("combining file %d %s" % (i,file))
                # f = nc.Dataset(os.path.join(ncdir, file))
                f = nc.Dataset(file)
                for j in range(gru_vars_num):
                    gru_var_name = gru_vars[j][0]
                    dim_index = gru_vars[j][1]
                    data=f[gru_var_name][:]
                    if i == 0:
                        Dict[gru_var_name]=data
                    else:
                        Dict[gru_var_name]=np.concatenate((Dict[gru_var_name],data),axis=dim_index)

                for j in range(hru_vars_num):
                    hru_var_name = hru_vars[j][0]
                    dim_index = hru_vars[j][1]
                    data=f[hru_var_name][:]
                    if i == 0:
                        Dict[hru_var_name]=data
                    else:
                        Dict[hru_var_name]=np.concatenate((Dict[hru_var_name],data),axis=dim_index)

            # assign values for gru and hru dimensioned variables
            for j in range(gru_vars_num):
                dst.variables[gru_vars[j][0]][:] = Dict[gru_vars[j][0]]
            for j in range(hru_vars_num):
                dst.variables[hru_vars[j][0]][:] = Dict[hru_vars[j][0]]

            #if missing HRUs, this is slow or broken
            if missing:
                new_index = np.append(dst["gru"].values,missgru)
                dst.reindex({"gru": new_index})
                dst.sel(gru=missgru)["gruId"] = missgru

                new_index = np.append(dst["hru"].values,misshru)
                dst.reindex({"hru": new_index})
                dst.sel(gru=misshru)["hruId"] = misshru

            # Remove grus and hrus that are not in the range
            cut_start = start_gru - start_gru_old
            cut_end = end_gru_old - end_gru
            dst.variables['gru'][:] = dst.variables['gru'][cut_start:]
            dst.variables['gru'][:] = dst.variables['gru'][:-cut_end]
            dst.variables['hru'][:] = dst.variables['hru'][cut_start:]
            dst.variables['hru'][:] = dst.variables['hru'][:-cut_end]

            # Temporarily create gruId from hruId
            #if gru_num == hru_num:
            #    gruId = dst.createVariable('gruId', dst['hruId'].datatype, ('gru',))
            #    gruId.long_name = "ID of group of response unit (GRU)"
            #    gruId.units = dst['hruId'].units
            #    dst.variables['gruId'][:] = dst.variables['hruId'][:]
            #else:
            #    print('Warning: gruId variable cannot be created since it has different size from hruId')

        print("wrote output: %s" % (ctdir + '/' + out_name))

    return out_name

# -- end functions

def process_group(g):
    size_g = int(np.ceil(gru_num / num_groups))
    start_gru = g * size_g + 1
    end_gru = min((g + 1) * size_g, gru_num)
    concatenate_files_in_range(outfilelist0, ctdir, start_gru, end_gru)

# -- end functions

if __name__ == "__main__":
    with mp.Pool(processes=mp.cpu_count()) as pool:
        pool.map(process_group, range(num_groups))