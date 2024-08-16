#!/usr/bin/env python
''' Create a vector cold state file for SUMMA from constant values'''
#
# Author:    Andy Wood, Feb 2017, modified by Ashley Van Beusekom, 2024
#
# Notes:  quick n dirty to generate constant initial states across a domain
#         all values hardwired, just gets HRU index from an existing parameter file
#         improvements:  could read another cold state file to get list of variables
#                        to populate; or read a metadata dictionary (names, types, etc)
#
#         no mapping required here -- but one could map another resolution vals, similar
#         to the param mapping scripts
#
#         check: nSoil and nSnow might have to be 'int' instead of 'int64' in output
#
# Requirements:  run with a python (eg miniconda) 2.7 that includes netCDF4
# Run as: python gen_coldstate.py <existing_attributeFile_with_hruId> <output_netCDF> <hru_type (int|int64|str)>
# =========================================================================

import sys
import os
import time
import getopt
import numpy as np
import netCDF4 as nc4
#import xarray as xr

testing = True
########################################################################
#                                Subroutines                           #
########################################################################

def usage():
    use = '''
    Usage: %s -[hv] [--help] [--verbose] <existing_inputfile_with_hruId> <output_netCDF> <hru_type (int|int64|str)>
            -h, --help      Show this help message and exit
            -v, --verbose   Enable verbose mode
    ''' % sys.argv[0]
    sys.stderr.write(use)
    sys.exit(1)

def getNetCDFData(fn, varname):
    """Read <varname> variables available to be mapped from NetCDF <fn> """
    f = nc4.Dataset(fn,'r')
    data = f.variables[varname][:]
    f.close()
#    ds = xr.open_dataset(fn)
#    data = ds[varname]
    return data

def getOutputPolyIDs(nc_file):
    outPolyIDs  = getNetCDFData(nc_file, 'hruId')    
    hru_elev = getNetCDFData(nc_file, 'elevation')
    hru_area = getNetCDFData(nc_file, 'HRUarea')
    gruIDs = getNetCDFData(nc_file, 'gruId')
    print("read output outPolyIds ('hruId') from example domain file")
    return outPolyIDs,hru_elev,hru_area, gruIDs


# write dom, hru, variables to netcdf output file
def writeNC_state_vars(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'dom','hru',),fill_value='-999.0')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write ngl, gru, variables to netcdf output file
def writeNC_state_varsGL(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'ngl','gru',),fill_value='-999.0')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn,  scalarv, midSoil, midToto, ifcToto, hrus, gruIDs, hru_type, ndom, ngl):    
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("writing output file")
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_hru = nc_out.createDimension('hru', len(hrus))
    dim_gru = nc_out.createDimension('gru', len(np.unique(gruIDs)))
    dim_scalarv = nc_out.createDimension('scalarv', scalarv)
    dim_midSoil = nc_out.createDimension('midSoil', midSoil)
    dim_midToto = nc_out.createDimension('midToto', midToto)
    dim_ifcToto = nc_out.createDimension('ifcToto', ifcToto)    
    dim_ndom = nc_out.createDimension('dom', ndom)
    dim_ngl = nc_out.createDimension('ngl', ngl)

    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20  # EC
        dim_str = nc_out.createDimension('strlen', max_strlen)
        hruId = nc_out.createVariable('hruId', 'S1', ('hru', 'strlen'),fill_value='-999')  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))     
        gruId = nc_out.createVariable('gruId', 'S1', ('gru', 'strlen'),fill_value='-999')
        gruId[:] = nc4.stringtochar(np.asarray(np.unique(gruIDs),
                                  dtype='S{}'.format(max_strlen)))

    elif hru_type == 'int64':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'i8', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')
        gruId = nc_out.createVariable('gruId', 'i8', ('gru', ),fill_value='-999')
        gruId[:] = np.unique(gruIDs)

    elif hru_type == 'int':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'int', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')
        gruId = nc_out.createVariable('gruId', 'int', ('gru', ),fill_value='-999')
        gruId[:] = np.unique(gruIDs)

    else:
        # not recognized
        sys.exit("ERROR, hru_type not recognized: must be str, int64, or int")

    # add attribute    
    hruId.long_name = 'USGS HUC12 ID'
    gruId.long_name = 'GRU ID'

    return nc_out
    # leave netcdf file open


import sys
import getopt

############################################
#                  Main                    #
############################################

if __name__ == '__main__':
    if testing:
        # hardwired for testing
        nc_example_name = 'attributes.nc'
        nc_out_name = 'coldstate_glac.nc'
        hru_type = 'int'

    else:
        try:
            (opts, args) = getopt.getopt(sys.argv[1:], 'hv', ['help', 'verbose'])
        except getopt.GetoptError as err:
            print(str(err))  # will print something like "option --f not recognized"
            usage()

        verbose = False
        for (opt, val) in opts:
            if opt in ('-h', '--help'):
                usage()
            elif opt in ('-v', '--verbose'):
                verbose = True
            else:
                print(f"Option {opt} not recognized")
                usage()

        if len(args) != 3:
            usage()
        nc_example_name = args[0]   # template file (param file, etc)
        nc_out_name = args[1]   # output cold-state file
        hru_type = args[2]   # 'int' or 'string'
 
    # hardwired to forcing formats (hru index rather than grid)
    outPolyIDs, hru_elev, hru_area, gruIDs = getOutputPolyIDs(nc_example_name)        
    nOutPolygonsHRU = len(outPolyIDs)
    nOutPolygonsGRU = len(gruIDs)

    # === now start to create the cold state variables using the variable template ===

    # settings 8 layer soil, 0 snow
    scalarv = 1
    midSoil = 8
    midToto = 8
    ifcToto = 9
    midIce = 3
    midLake = 0
    midSnow = 0
    dT = 3600 # timestep of forcings in seconds
    lyrDepth  = [0.025, 0.075, 0.15, 0.25, 0.5, 0.5, 1.0, 1.5]
    lyrHeight = [0.0, 0.025, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0]

    ndom = 3
    ngl = 2

    # initialize netcdf file by storing dimensions and hru variable
    nc_out = writeNC_dims(nc_out_name, scalarv, midSoil, midToto, ifcToto,
                        outPolyIDs, gruIDs, hru_type, ndom, ngl)

    # === now loop through variables and write
    #  this could be done by looping through the input state file and xferring values
    #  domain order has to be upland, glacier ablation, glacier accumulation, wetland
    #   has to have upland, then if has glacier has to have ablation and accumulation

    # nSoil, nSnow, nIce, nLake
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), midSoil, dtype='f8')
    writeNC_state_vars(nc_out, 'nSoil', 'scalarv', 'f8', newVarVals)
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), midSnow, dtype='f8')
    writeNC_state_vars(nc_out, 'nSnow', 'scalarv', 'f8', newVarVals)
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), midIce,  dtype='f8')
    writeNC_state_vars(nc_out, 'nIce', 'scalarv', 'f8', newVarVals)
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), midLake, dtype='f8')
    writeNC_state_vars(nc_out, 'nLake', 'scalarv', 'f8', newVarVals)

    # dT
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), dT, dtype='f8')
    writeNC_state_vars(nc_out, 'dt_init', 'scalarv', 'f8', newVarVals)

    # area
    newVarVals = np.full((1, ndom, nOutPolygonsHRU), hru_area / 3, dtype='f8')
    writeNC_state_vars(nc_out, 'DOMarea', 'scalarv', 'f8', newVarVals)

    # elevation
    newVarVals = np.full((1, nOutPolygonsHRU), hru_elev, dtype='f8')
    writeNC_state_vars(nc_out, 'DOMelev', 'scalarv', 'f8', newVarVals)

    # SWE, SnowDepth, SfcMeltPond, SnowAlbedo, CanopyLiq, CanopyIce
    newVarVals = np.full((1, ndom, nOutPolygonsGRU), 0.0, dtype='f8')
    writeNC_state_vars(nc_out, 'scalarSWE', 'scalarv', 'f8', newVarVals)
    writeNC_state_vars(nc_out, 'scalarSnowDepth', 'scalarv', 'f8', newVarVals)
    writeNC_state_vars(nc_out, 'scalarSfcMeltPond', 'scalarv', 'f8', newVarVals)
    writeNC_state_vars(nc_out, 'scalarSnowAlbedo', 'scalarv', 'f8', newVarVals)
    writeNC_state_vars(nc_out, 'scalarCanopyLiq', 'scalarv', 'f8', newVarVals)

    # glacier area
    newVarVals = np.full((1,ngl,nOutPolygonsHRU), 0.0, dtype='f8')
    writeNC_state_varsGL(nc_out, 'glacAblArea', 'scalarv', 'f8', newVarVals)     
    writeNC_state_varsGL(nc_out, 'glacAccArea', 'scalarv', 'f8', newVarVals)

    # CanairTemp, CanopyTemp
    newVarVals = np.full((1,ndom, nOutPolygonsHRU), 283.16, dtype='f8')     
    writeNC_state_vars(nc_out, 'scalarCanairTemp', 'scalarv', 'f8', newVarVals)
    writeNC_state_vars(nc_out, 'scalarCanopyTemp', 'scalarv', 'f8', newVarVals)

    # AquiferStorage
    newVarVals = np.full((1,ndom,nOutPolygonsHRU), 1.0, dtype='f8')      
    writeNC_state_vars(nc_out, 'scalarAquiferStorage', 'scalarv', 'f8', newVarVals)

    # layer MatricHead
    newVarVals = np.full((midSoil,ndom,nOutPolygonsHRU), -1.0, dtype='f8')
    writeNC_state_vars(nc_out, 'mLayerMatricHead', 'midSoil', 'f8', newVarVals)

    # layer Temp
    newVarVals = np.full((midToto,ndom,nOutPolygonsHRU), 283.16, dtype='f8')
    writeNC_state_vars(nc_out, 'mLayerTemp', 'midToto', 'f8', newVarVals)

    # layer VolFracLiq
    newVarVals = np.full((midToto,ndom,nOutPolygonsHRU), 0.2, dtype='f8')
    writeNC_state_vars(nc_out, 'mLayerVolFracLiq', 'midToto', 'f8', newVarVals)

    # layer VolFracIce
    newVarVals = np.full((midToto,ndom,nOutPolygonsHRU), 0.0, dtype='f8')
    writeNC_state_vars(nc_out, 'mLayerVolFracIce', 'midToto', 'f8', newVarVals)

    # layer Depth, Height
    newVarVals = np.full((nOutPolygonsHRU,ndom,midToto), lyrDepth, dtype='f8').transpose()
    writeNC_state_vars(nc_out, 'mLayerDepth', 'midToto', 'f8', newVarVals)
    newVarVals = np.full((nOutPolygonsHRU,ndom,ifcToto), lyrHeight, dtype='f8').transpose()
    writeNC_state_vars(nc_out, 'iLayerHeight', 'ifcToto', 'f8', newVarVals)        

    nc_out.close()


