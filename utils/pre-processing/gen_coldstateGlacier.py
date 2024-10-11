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
from netCDF4 import Dataset
import shutil
#import xarray as xr

testing = True
glac_dom = True 
wtld_dom = False
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
    hru2gru = getNetCDFData(nc_file, 'hru2gruId')
    print("read data from attribute file")
    return outPolyIDs,hru_elev,hru_area, gruIDs,hru2gru

# write gru, variables to netcdf output file
def writeNC_state_vars_GRU(nc_out, newVarName, newVarType, newVarVals):
    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("adding attribute data")
    ncvar = nc_out.createVariable(newVarName, newVarType, ('gru',),fill_value='-999')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write hru, variables to netcdf output file
def writeNC_state_vars_HRU(nc_out, newVarName, newVarType, newVarVals):
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("adding attribute data")
    ncvar = nc_out.createVariable(newVarName, newVarType, ('hru',),fill_value='-999')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write dom, hru, variables to netcdf output file
def writeNC_state_vars_HRU_DOM(nc_out, newVarName, newVarDim, newVarType, newVarVals):
    """ Write <vars>[hru dom] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("adding HRU_DOM data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'hru','dom',),fill_value='-999.0')   
    ncvar[:] = newVarVals   # store data in netcdf file

# write ngl, gru, variables to netcdf output file
def writeNC_state_vars_GRU_NGL(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding GRU_NGL data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'gru',),fill_value='-999.0')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn,  scalarv, midSoil, midToto, ifcToto, hrus, grus, hru_type, ndom, ngl):    
    """ Write <vars>[hru] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("writing output file")
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_hru = nc_out.createDimension('hru', len(hrus))
    dim_gru = nc_out.createDimension('gru', len(grus))
    dim_scalarv = nc_out.createDimension('scalarv', scalarv)
    dim_midSoil = nc_out.createDimension('midSoil', midSoil)
    dim_midToto = nc_out.createDimension('midToto', midToto)
    dim_ifcToto = nc_out.createDimension('ifcToto', ifcToto)    
    dim_ndom = nc_out.createDimension('dom', ndom) # max number of domains in any HRU
    dim_ngl = nc_out.createDimension('ngl', ngl) # max number of glaciers in any GRU

    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20  # EC
        dim_str = nc_out.createDimension('strlen', max_strlen)
        hruId = nc_out.createVariable('hruId', 'S1', ('hru', 'strlen'),fill_value='-999')  
        hruId[:] = nc4.stringtochar(np.asarray(hrus,
                                  dtype='S{}'.format(max_strlen)))     
        gruId = nc_out.createVariable('gruId', 'S1', ('gru', 'strlen'),fill_value='-999')
        gruId[:] = nc4.stringtochar(np.asarray(np.unique(grus),
                                  dtype='S{}'.format(max_strlen)))

    elif hru_type == 'int64':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'i8', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')
        gruId = nc_out.createVariable('gruId', 'i8', ('gru', ),fill_value='-999')
        gruId[:] = grus

    elif hru_type == 'int':
        # integer HRU
        hruId = nc_out.createVariable('hruId', 'int', ('hru', ),fill_value='-999')   
        hruId[:] = hrus
        #hruId[:] = np.asarray(hrus, dtype='int')
        gruId = nc_out.createVariable('gruId', 'int', ('gru', ),fill_value='-999')
        gruId[:] = grus

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
        nc_attribute_name = 'attributes.nc'
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
        nc_attribute_name = args[0]     # attribute file with HRU index
        nc_out_name = args[1]           # output cold-state file
        hru_type = args[2]              # 'int' or 'string'
 
    # hardwired to forcing formats (hru index rather than grid)
    outPolyIDs, hru_elev, hru_area, gruIDs, hru2gru = getOutputPolyIDs(nc_attribute_name)        
    nOutPolygonsHRU = len(outPolyIDs)
    nOutPolygonsGRU = len(gruIDs)

    # === create the modified attribute netcdf file ===
    nc_attribute_name_new = nc_attribute_name[:-3] + '_glac.nc'
    shutil.copy(nc_attribute_name, nc_attribute_name_new)
    print("created new attribute file")

    # add nGlacier and nWetland to attribute file
    nGlacier0 = 0
    nWetland0 = 0
    if glac_dom: # add glaciers to every GRU for testing, SHOULD BE READ FROM FILE
        nGlacier0 = 2 
    if wtld_dom: # add wetlands to every HRU for testing, SHOULD BE READ FROM FILE
        nWetland0 = 1

    with Dataset(nc_attribute_name_new, 'a') as nc_out:
        nGlacier = np.full((nOutPolygonsGRU), nGlacier0, dtype='i4')
        writeNC_state_vars_GRU(nc_out, 'nGlacier', 'i4', nGlacier)
        nWetland = np.full((nOutPolygonsGRU), nWetland0, dtype='i4')
        writeNC_state_vars_GRU(nc_out, 'nWetland', 'i4', nWetland)

    ngl = max(nGlacier)

    # === now start to create the cold state variables using the variable template ===
    dT = 3600 # timestep of forcings in seconds

    # settings 8 layer soil, 0 snow
    scalarv = 1
    midSoil = 8
    midSoil_glac = 0
    midSoil_wtld = 0
    midLake = 0
    midIce = 0
    midSnow = 0
    ndom = 1
    midToto = 8
    midToto_glac = 0
    midToto_wtld = 0
    lyrHeight0 = [0.0, 0.025, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0]
    lyrHeight0_np = np.array(lyrHeight0)
    lyrDepth0 = lyrHeight0_np[1:] - lyrHeight0_np[:-1]

    # domain order here is 1)upland, 2)glacier accumulation, 3)glacier ablation, 4)wetland
    #   adjust these values as needed SHOULD READ FROM A FILE
    upld_frac = 1.0
    glac_frac = 0.5
    wtld_frac = 0.25
    indxWtld = 1
    if glac_dom: 
        ndom += 2
        midIce = 3
        midSoil_glac = 2
        midToto_glac = 7
        lyrHeight_glac = [0.0, 0.025, 0.125, 0.275, 0.575, 2.375, 7.125, 30.125]
        lyrHeight_glacnp = np.array(lyrHeight_glac)
        lyrDepth_glac = lyrHeight_glacnp[1:] - lyrHeight_glacnp[:-1]
        indxWtld += 2
        upld_frac = upld_frac - glac_frac
    if wtld_dom: 
        ndom += 1
        midLake = 5
        midSoil_wtld = 5
        midToto_wtld = 9
        lyrHeight_wtld = [0.0, 0.05, 0.2, 0.5, 1, 2.0, 3.0, 4.0, 5.0, 8.0]
        lyrHeight_wtldnp = np.array(lyrHeight_wtld)
        lyrDepth_wtld = lyrHeight_wtldnp[1:] - lyrHeight_wtldnp[:-1]
        upld_frac = upld_frac - wtld_frac

    midToto = max(midToto_glac, midToto_wtld,midToto)
    midSoil = max(midSoil_glac, midSoil_wtld,midSoil)
    ifcToto = midToto + 1

    # initialize layer variables
    toto0       = np.full((midToto, nOutPolygonsHRU, ndom), 0.0, dtype='f8')
    totopoint2  = np.full((midToto, nOutPolygonsHRU, ndom), 0.2, dtype='f8')
    soilneg1    = np.full((midSoil, nOutPolygonsHRU, ndom), -1.0, dtype='f8')
    toto283     = np.full((midToto, nOutPolygonsHRU, ndom), 283.16, dtype='f8')

    scalar0     = np.full((1, nOutPolygonsHRU, ndom), 0.0, dtype='f8')
    scalar283   = np.full((1, nOutPolygonsHRU, ndom), 283.16, dtype='f8')
    scalar1     = np.full((1, nOutPolygonsHRU, ndom), 1.0, dtype='f8')

    # layer Depth, Height, layer types adjust for glacier and wetland
    lyrDepth = np.zeros((ndom, nOutPolygonsHRU, midToto), dtype='f8')
    lyrDepth[0,:,0:len(lyrDepth0)]  = lyrDepth0
    lyrHeight = np.zeros((ndom, nOutPolygonsHRU, ifcToto), dtype='f8')
    lyrHeight[0,:,0:len(lyrHeight0)] = lyrHeight0
    midSoil_dom = np.full((1, nOutPolygonsHRU, ndom), midSoil, dtype='f8')
    midIce_dom = scalar0.copy()
    midLake_dom = scalar0.copy()
    dom_area = np.full((1, nOutPolygonsHRU, ndom), hru_area, dtype='f8')
    dom_area[0,0,:] = hru_area * upld_frac
    dom_elev = np.full((1, nOutPolygonsHRU, ndom), hru_elev, dtype='f8')

    domType = np.full((1, nOutPolygonsHRU, ndom), 1, dtype='i4')
    if glac_dom: # NOTE, if HRU glacier area is 0, midIce_dom should be 0
        lyrDepth[1,:,0:len(lyrDepth_glac)] = lyrDepth_glac
        lyrDepth[2,:,0:len(lyrDepth_glac)] = lyrDepth_glac
        lyrHeight[1,:,0:len(lyrHeight_glac)] = lyrHeight_glac
        lyrHeight[2,:,0:len(lyrHeight_glac)] = lyrHeight_glac
        midIce_dom[:,0,1] = midIce
        midIce_dom[:,0,2] = midIce
        midSoil_dom[:,0,1] = midSoil_glac
        midSoil_dom[:,0,2] = midSoil_glac
        domType[:,0,1] = 2 # glacier accumulation
        dom_area[:,0,1] = hru_area * glac_frac*0.5  # glacier accumulation 50% of glacier, SHOULD BE READ FROM FILE
        dom_elev[:,0,1] = hru_elev * 1.5            # for testing, SHOULD BE READ FROM FILE
        domType[:,0,2] = 3 # glacier ablation
        dom_area[:,0,2] = hru_area * glac_frac*0.5  # glacier ablation 50% of glacier, SHOULD BE READ FROM FILE
        dom_elev[:,0,2] = hru_elev * 0.5            # for testing, SHOULD BE READ FROM FILE

    if wtld_dom: # NOTE, if HRU wetland area is 0, midLake_dom should be 0
        lyrDepth[indxWtld,:,0:len(lyrDepth_wtld)] = lyrDepth_wtld
        lyrHeight[indxWtld,:,0:len(lyrHeight_wtld)] = lyrHeight_wtld
        midLake_dom[:,0,indxWtld] = midLake
        midSoil_dom[:,0,indxWtld] = midSoil_wtld
        domType[:,0,indxWtld] = 4
        dom_area[:,0,indxWtld] = hru_area * wtld_frac 
        dom_elev[:,0,indxWtld] = hru_elev           # assume wetland elev same as upland

    lyrDepth = lyrDepth.transpose()
    lyrHeight = lyrHeight.transpose()

    # initialize netcdf file by storing dimensions and hru variable
    nc_out = writeNC_dims(nc_out_name, scalarv, midSoil, midToto, ifcToto,
                        outPolyIDs, gruIDs, hru_type, ndom, ngl)

    # === now loop through variables and write
    #  this could be done by looping through the input state file and xferring values
    #  domain order has to be upland, glacier ablation, glacier accumulation, wetland
    #   has to have upland, then if has glacier has to have ablation and accumulation

    # domType
    writeNC_state_vars_HRU_DOM(nc_out, 'domType', 'scalarv', 'f8', domType)

    # layer Depth, Height
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerDepth', 'midToto', 'f8', lyrDepth)        # Depth
    writeNC_state_vars_HRU_DOM(nc_out, 'iLayerHeight', 'ifcToto', 'f8', lyrHeight)      # Height

    # nSoil, nSnow, nIce, nLake
    writeNC_state_vars_HRU_DOM(nc_out, 'nSoil', 'scalarv', 'f8', midSoil_dom)           # nSoil
    writeNC_state_vars_HRU_DOM(nc_out, 'nSnow', 'scalarv', 'f8', scalar0 )              # nSnow start at 0
    writeNC_state_vars_HRU_DOM(nc_out, 'nIce', 'scalarv', 'f8', midIce_dom)             # nIce
    writeNC_state_vars_HRU_DOM(nc_out, 'nLake', 'scalarv', 'f8', midLake_dom)           # nLake

    # dT
    newVarVals = np.full((ndom, 1, nOutPolygonsHRU), dT, dtype='f8')
    writeNC_state_vars_HRU_DOM(nc_out, 'dt_init', 'scalarv', 'f8', newVarVals)

    # area and elevation
    writeNC_state_vars_HRU_DOM(nc_out, 'DOMarea', 'scalarv', 'f8', dom_area)           # DOMarea
    writeNC_state_vars_HRU_DOM(nc_out, 'DOMelev', 'scalarv', 'f8', dom_elev)           # DOMelev

    # SWE, SnowDepth, SfcMeltPond, SnowAlbedo, CanopyLiq, CanopyIce
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSWE', 'scalarv', 'f8', scalar0)           # SWE
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSnowDepth', 'scalarv', 'f8', scalar0)     # SnowDepth
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSfcMeltPond', 'scalarv', 'f8', scalar0)   # SfcMeltPond
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSnowAlbedo', 'scalarv', 'f8', scalar0)    # SnowAlbedo
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyLiq', 'scalarv', 'f8', scalar0)     # CanopyLiq
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyIce', 'scalarv', 'f8', scalar0)     # CanopyIce

    # glacier area, just divide by ngl for testing, SHOULD BE READ FROM FILE
    totAccArea = np.zeros(nOutPolygonsGRU)
    totAblArea = np.zeros(nOutPolygonsGRU)
    accArea = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    ablArea = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    for i,g in enumerate(gruIDs):
        totAccArea[i] = dom_area[0,hru2gru==g,1].sum()
        totAblArea[i] = dom_area[0,hru2gru==g,2].sum()
        if nGlacier[i]>0: 
            accArea[0,:nGlacier[i],i] = totAccArea[i]/nGlacier[i]
            ablArea[0,:nGlacier[i],i] = totAblArea[i]/nGlacier[i]
    

    newVarVals = np.full((ngl,nOutPolygonsGRU), ablArea, dtype='f8')
    writeNC_state_vars_GRU_NGL(nc_out, 'glacAblArea', 'ngl', 'f8', newVarVals)
    newVarVals = np.full((ngl,nOutPolygonsGRU), accArea, dtype='f8')     
    writeNC_state_vars_GRU_NGL(nc_out, 'glacAccArea', 'ngl', 'f8', newVarVals)

    # scalar CanairTemp, CanopyTemp, AquiferStorage
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanairTemp', 'scalarv', 'f8', scalar283)  # CanairTemp, does not exist for glacier or lake
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyTemp', 'scalarv', 'f8', scalar283)  # CanopyTemp, does not exist for glacier or lake
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarAquiferStorage', 'scalarv', 'f8',scalar1) # AquiferStorage, does not exist for glacier

    # layer MatricHead, Temp, VolFracLiq, VolFracIce
    if glac_dom:
        for i in range(ndom): # put ice layers at -5 C and all ice similar to Giese et al. 2020, otherwise need to spin up 40 yrs
            if i>0 and i<3: 
                toto283[   midSoil_glac:(midIce+midSoil_glac),:,i] = 268.16 # or 273.16?
                toto0[     midSoil_glac  ,:,i] = 0.90 # could be 0.9 per Bradford et al. 2009, less air as go deeper
                toto0[     midSoil_glac+1,:,i] = 0.91
                toto0[     midSoil_glac+2,:,i] = 0.93
                toto0[     midSoil_glac+3,:,i] = 0.95 
                toto0[     midSoil_glac+4,:,i] = 0.98
                totopoint2[midSoil_glac:(midIce+midSoil_glac),:,i] = 0.0
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerMatricHead', 'midSoil', 'f8', soilneg1)   # MatricHead
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerTemp', 'midToto', 'f8', toto283)          # Temp
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerVolFracLiq', 'midToto', 'f8', totopoint2) # VolFracLiq
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerVolFracIce', 'midToto', 'f8', toto0)      # VolFracIce


    nc_out.close()


