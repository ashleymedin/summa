#!/usr/bin/env python
''' Create a vector cold state file for SUMMA glacier or wetland from constant values'''
#
# Author: Ashley Van Beusekom, 2024
#
# Notes:  quick n dirty to generate constant initial states across a domain
#         all values hardwired, just gets HRU index from an existing parameter file
#
#         no mapping required here -- but one could map another resolution vals, similar
#         to the param mapping scripts
##
# Requirements:  run with a python (eg miniconda) 2.7 that includes netCDF4
# NEED TO RUN AFTER gen_glacTopo.py
# Run as: python gen_coldstateGlacier.py <existing_attributeFile_with_hruId> <output_netCDF> <hru_type (int|int64|str)>
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
thick4area = 0.1 #thickness of ice to count as glacier area
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

def getOutputGlac_bed(nc_file):
    nGlacier  = getNetCDFData(nc_file, 'nGlacier')
    gruIDs = getNetCDFData(nc_file, 'gruId')
    gridIDs  = getNetCDFData(nc_file, 'gridId')
    nx = getNetCDFData(nc_file, 'nx')
    ny = getNetCDFData(nc_file, 'ny')
    dx  = getNetCDFData(nc_file, 'dx')
    dy  = getNetCDFData(nc_file, 'dy')
    cell2hruId  = getNetCDFData(nc_file, 'cell2hruId')
    bed_elev  = getNetCDFData(nc_file, 'bed_elev')
    debris_thick  = getNetCDFData(nc_file, 'debris_thick')
    stage = getNetCDFData(nc_file, 'stage')
    print("read data from glacier bed file")
    return nGlacier,gruIDs,gridIDs,nx,ny,dx,dy,cell2hruId,bed_elev,debris_thick,stage

def getOutputGlac_surf(nc_file):
    gruIDs = getNetCDFData(nc_file, 'gruId')
    gridIDs  = getNetCDFData(nc_file, 'gridId')
    AAR = getNetCDFData(nc_file, 'AAR')
    surface_elev  = getNetCDFData(nc_file, 'surface_elev')
    print("read data from glacier surface file")
    return gruIDs,gridIDs,AAR,surface_elev

# write gru variables to netcdf output file
def writeNC_state_vars_GRU(nc_out, newVarName, newVarType, newVarVals):
    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("adding attribute data")
    ncvar = nc_out.createVariable(newVarName, newVarType, ('gru',),fill_value='-999')    
    ncvar[:] = newVarVals   # store data in netcdf file

# write dom, hru, variables to netcdf output file
def writeNC_state_vars_HRU_DOM(nc_out, newVarName, newVarDim, newVarType, newVarVals):
    """ Write <vars>[hru dom] array in netCDF4 file,<fn> and variable of
        <varname> """
    print("adding HRU_DOM data")
    ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'hru','dom',),fill_value='-999.0')   
    ncvar[:] = newVarVals   # store data in netcdf file

# write gru nonscalar, variables to netcdf output file
def writeNC_state_vars_GRU_VEC(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding GRU_VEC data")
    if newVarType=='i4' or newVarType=='i8':
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'gru',),fill_value='-999')  
    else:
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'gru',),fill_value='-999.0')  
    ncvar[:] = newVarVals   # store data in netcdf file

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn,  scalarv, midSoil, midToto, ifcToto, hrus, grus, hru_type, ndom, glac):    
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
    dim_ngl = nc_out.createDimension('glac', glac) # max number of glaciers in any GRU

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

############################################
#                  Main                    #
############################################

if __name__ == '__main__':
    if testing:
        # hardwired for testing
        nc_attribute_name = 'attributes.nc'
        nc_out_name0 = 'coldstate.nc'
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
        nc_out_name0 = args[1]           # output cold-state file
        hru_type = args[2]              # 'int' or 'string'
 
    # hardwired to forcing formats (hru index rather than grid)
    outPolyIDs, hru_elev, hru_area, gruIDs, hru2gru = getOutputPolyIDs(nc_attribute_name)        
    nOutPolygonsHRU = len(outPolyIDs)
    nOutPolygonsGRU = len(gruIDs)
    nc_out_name = nc_out_name0[:-3] + '_glac.nc'

    # === read glacier data ===
    if glac_dom:
        nc_glacier_name = nc_attribute_name[:-3] + '_glacBedTopo.nc'
        nc_glacierS_name = nc_out_name0[:-3] + '_glacSurfTopo.nc'
        nGlacier, gruIDs2, gridIDs, nx, ny, dx, dy, cell2hruId, bed_elev, debris_thick0, stage = getOutputGlac_bed(nc_glacier_name)
        gruIDs3, gridIDs2, AAR, surface_elev = getOutputGlac_surf(nc_glacierS_name)
        glac = max(nGlacier)

        # Mapping to gruIDs
        gruID_map = {gru: i for i, gru in enumerate(gruIDs)}

        # Remap gruIDs2 and gruIDs3 to match gruIDs
        gruIDs2_mapped = np.array([gruID_map[gru] for gru in gruIDs2])
        gruIDs3_mapped = np.array([gruID_map[gru] for gru in gruIDs3])

        # Remap all variables to match gruIDs
        nGlacier = nGlacier[gruIDs2_mapped]
        gridIDs  = gridIDs[:,gruIDs2_mapped]
        nx = nx[:,gruIDs2_mapped]
        ny = ny[:,gruIDs2_mapped]
        dx = dx[:,gruIDs2_mapped]
        dy = dy[:,gruIDs2_mapped]
        cell2hruId   = cell2hruId[:,:,:,gruIDs2_mapped]
        bed_elev     = bed_elev[:,:,:,gruIDs2_mapped]
        debris_thick0 = debris_thick0[:,gruIDs2_mapped]
        stage        = stage[:,gruIDs2_mapped]
        gridIDs2     = gridIDs2[:,gruIDs3_mapped]
        AAR          = AAR[:,gruIDs3_mapped]
        surface_elev = surface_elev[:,:,:,gruIDs3_mapped]

        # Remap gridIDs2 to gridIDs for each gru
        for i, g in enumerate(gruIDs):
            gridIDs_map = {glac: j for j, glac in enumerate(gridIDs[:,i])}
            gridIDs2_mapped = np.array([gridIDs_map[glac] for glac in gridIDs2[:,i]])

            # Remap all variables to match gridIDs
            AAR[:,i] = AAR[gridIDs2_mapped,i]
            surface_elev[:,:,:,i] = surface_elev[:,:,gridIDs2_mapped,i]

    else:
        nGlacier = np.zeros(nOutPolygonsGRU, dtype='i4')

    # === create the modified attribute netcdf file ===
    nc_attribute_name_new = nc_attribute_name[:-3] + '_glac.nc'
    shutil.copy(nc_attribute_name, nc_attribute_name_new)
    print("created new attribute file")

    # add wetlands to every HRU for testing, SHOULD BE READ FROM FILE
    nWetland0 = 0
    if wtld_dom: 
        nWetland0 = 1

    with Dataset(nc_attribute_name_new, 'a') as nc_out:
        writeNC_state_vars_GRU(nc_out, 'nGlacier', 'i4', nGlacier)
        nWetland = np.full((nOutPolygonsGRU), nWetland0, dtype='i4')
        writeNC_state_vars_GRU(nc_out, 'nWetland', 'i4', nWetland)

    # === now start to create the cold state variables using the variable template ===
    dT = 3600 # timestep of forcings in seconds

    # settings 8 layer soil, 0 snow
    scalarv = 1
    midSoil = 8
    midSoil_glac = np.zeros((nOutPolygonsHRU,2), dtype='i4') # [no debris, debris]
    midSoil_wtld = np.zeros(nOutPolygonsHRU,dtype='i4')
    midLake = 0     # assume same number of ice layers for every glacier
    midGlce = 0      # assume same number of lake layers for every lake
    midSnow = 0
    ndom = 1
    midToto = 8
    midToto_glac = np.zeros((nOutPolygonsHRU,2), dtype='i4') # [no debris, debris]
    midToto_wtld = np.zeros(nOutPolygonsHRU,dtype='i4')
    lyrHeight0 = [0.0, 0.025, 0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 4.0]
    lyrHeight0_np = np.array(lyrHeight0)
    lyrDepth0 = lyrHeight0_np[1:] - lyrHeight0_np[:-1]

    # domain order here is 1)upland, 2)glacier accumulation, 3)glacier ablation clean, 4)glacier ablation debris, 5)wetland
    #   adjust these values as needed SHOULD READ FROM A FILE
    upld_frac = np.ones(nOutPolygonsHRU, dtype='f8')
    glacAcc_frac = np.zeros(nOutPolygonsHRU,dtype='f8')
    glacCln_frac = np.zeros(nOutPolygonsHRU, dtype='f8')
    wtld_frac = np.zeros(nOutPolygonsHRU, dtype='f8')
    upld_elev = np.full(nOutPolygonsHRU,hru_elev, dtype='f8')
    glacAcc_elev = np.zeros(nOutPolygonsHRU, dtype='f8')
    glacCln_elev = np.zeros(nOutPolygonsHRU, dtype='f8')
    wtld_elev = np.zeros(nOutPolygonsHRU, dtype='f8')
    debris_thick = np.zeros(nOutPolygonsHRU, dtype='f8')
    glacDbr_frac = np.zeros(nOutPolygonsHRU, dtype='f8')
    glacDbr_elev = np.zeros(nOutPolygonsHRU, dtype='f8')
    glacNoDeb_frac = np.zeros(nOutPolygonsHRU, dtype='f8')
    glacNoDeb_elev = np.zeros(nOutPolygonsHRU, dtype='f8')
    indWtld = 1 # initialize, will be incremented for other domains
    debr_dom = 0 # initialize to no debris
    clean_dom = 0 # initialize to no clean ice
    if glac_dom: 
        ablArea = np.zeros((glac,nOutPolygonsGRU), dtype='f8')
        accArea = np.zeros((glac,nOutPolygonsGRU), dtype='f8')
        totVolume = np.zeros(nOutPolygonsGRU, dtype='f8')
        for i,g in enumerate(gruIDs):
            for j in range(nGlacier[i]):
                bed_elev0 = bed_elev[:,:,j,i]
                surface_elev0 = surface_elev[:,:,j,i]
                hgt = surface_elev0 - bed_elev0
                cell2hruId0 = cell2hruId[:,:,j,i]
                dly = dy[j,i]*dx[j,i]
                
                # calculate area and volume
                glac_area = np.sum(np.where((hgt>thick4area), dly, 0))
                accArea[j,i] = glac_area * AAR[j,i]
                ablArea[j,i] = glac_area - accArea[j,i]
                glac_vol = np.sum(np.where((hgt>0), hgt*dly, 0))    # m^9 of ice
                totVolume[i] = totVolume[i] + glac_vol

                # calculate the ELA
                # Sort the elevations and corresponding areas from lowest elev to highest
                sorted_indices = np.argsort(surface_elev0[hgt>0].flatten())
                sorted_elevations = surface_elev0[hgt>0].flatten()[sorted_indices]
                # Initialize the accumulated area
                accumulated_area = 0
                ELA_elev = -999.0
                # Loop through the sorted elevations and areas
                for idy in range(len(sorted_elevations)):
                    accumulated_area += dly
                    if accumulated_area >= ablArea[j, i]:
                        ELA_elev = sorted_elevations[idy]
                        print('Glacier ', gridIDs[j, i], ' ELA, area, vol:', ELA_elev, glac_area, glac_vol)
                        break
                if ELA_elev == -999.:
                    print('Glacier ', gridIDs[j,i],' ELA not found')
                
                # calculate area and elevation by HRU
                for k,h in enumerate(outPolyIDs):
                    if hru2gru[k] == g:
                        debris_thick[k] = debris_thick[k] + np.sum(np.where((cell2hruId0==h) & (hgt>0) & (surface_elev0 <ELA_elev), debris_thick0[j,i]*dly*stage[j,i], 0)) # debris thickness by glacier in GRU, average by HRU
                        glacCln_frac[k] = glacCln_frac[k] + np.sum(np.where((cell2hruId0==h) & (hgt>0) & (surface_elev0 <ELA_elev), dly*(1.0-stage[j,i]), 0))/hru_area[k]   
                        glacDbr_frac[k] = glacDbr_frac[k] + np.sum(np.where((cell2hruId0==h) & (hgt>0) & (surface_elev0 <ELA_elev), dly*stage[j,i], 0))/hru_area[k] 
                        glacAcc_frac[k] = glacAcc_frac[k] + np.sum(np.where((cell2hruId0==h) & (hgt>0) & (surface_elev0>=ELA_elev), dly, 0))/hru_area[k]
                        glacCln_elev[k] = glacCln_elev[k] + surface_elev0[(cell2hruId0==h) & (hgt>0) & (surface_elev0< ELA_elev)].mean()/nGlacier[i]
                        glacDbr_elev[k] = glacCln_elev[k] # make all ablation elevations the same for now
                        glacAcc_elev[k] = glacAcc_elev[k] + surface_elev0[(cell2hruId0==h) & (hgt>0) & (surface_elev0>=ELA_elev)].mean()/nGlacier[i]
                        glacDbr_frac[k] = glacDbr_frac[k] + np.sum(np.where((cell2hruId0==h) & (hgt>0) & (surface_elev0 <ELA_elev), dly*stage[j,i], 0))/hru_area[k] 

        ndom += 2
        midGlce = 5
        midSoil_glac0 = 2 # if has debris, will have this many soil layers
        midToto_glac0 = midGlce + midSoil_glac0 # if has debris, will have this many total layers
        lyrHeight_glac = np.zeros((nOutPolygonsHRU,midToto_glac0+1,2), dtype='f8') # [no debris, debris]
        lyrHeight_glacNoDeb = np.zeros((nOutPolygonsHRU,midToto_glac0+1,2), dtype='f8')
        for k,h in enumerate(outPolyIDs):
            if glacDbr_frac[k]>0.0:
                debris_thick[k] = debris_thick[k]/(glacDbr_frac[k]*hru_area[k])
                midSoil_glac[k,1] = midSoil_glac0 # midSoil_glac[k,0] always =0
                for layer in range(1, midSoil_glac[k,1] + 1):
                    lyrHeight_glac[k, layer,1] = debris_thick[k] * (layer ** 2 / midSoil_glac[k,1] ** 2)
                ice_depth = 30.0
                ice_layDepth = [0.15, 0.45, 2.25, 7.0, 30.0]
                lyrHeight_glac[k, 1:midGlce+1,0] = ice_layDepth
                for layer in range(midSoil_glac[k,1]+1, midGlce+midSoil_glac[k,1]+1): # midSoil_glac[k,1] = 0 if no debris in HRU, or if no debri domain
                    lyrHeight_glac[k, layer,1] = debris_thick[k] + ice_layDepth[layer-midSoil_glac[k,1]-1]

            midToto_glac[k,:] = midGlce+midSoil_glac[k,:]
                
        lyrHeight_glacnp = np.array(lyrHeight_glac)
        lyrDepth_glac = lyrHeight_glacnp[:,1:,:] - lyrHeight_glacnp[:,:-1,:]
        indWtld += 2
        if np.sum(stage) > 0.0: # if there is debris
            if np.sum(stage)<nGlacier: 
                ndom += 1 #
                indDebr = 3
                indWtld += 1
                clean_dom = 1
            else: # no clean ice on any glacier
                indDebr = 2
            debr_dom = 1

    if wtld_dom: 
        ndom += 1
        # These should be by HRU
        midLake = 4
        midSoil_wtld = 5
        midToto_wtld = midLake + midSoil_wtld
        lyrHeight_wtld = [0.0, 0.05, 0.2, 0.5, 1, 2.0, 3.0, 4.0, 5.0, 8.0]
        lyrHeight_wtldnp = np.array(lyrHeight_wtld)
        lyrDepth_wtld = lyrHeight_wtldnp[1:] - lyrHeight_wtldnp[:-1]
        wtld_frac = 0.25 # SHOULD BE READ FROM FILE
        wtld_elev = hru_elev # SHOULD BE READ FROM FILE

    upld_frac = upld_frac - glacCln_frac - glacDbr_frac - glacAcc_frac - wtld_frac
    upld_elev = upld_elev - glacCln_frac * glacCln_elev - glacDbr_frac * glacDbr_elev - glacAcc_frac * glacAcc_elev - wtld_frac * wtld_elev
    upld_elev = upld_elev/upld_frac

    midToto = max(np.max(midToto_glac), midToto_wtld,midToto)
    midSoil = max(np.max(midSoil_glac), midSoil_wtld,midSoil)
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
    midGlce_dom = scalar0.copy()
    midLake_dom = scalar0.copy()
    dom_area = np.full((1, nOutPolygonsHRU, ndom), hru_area, dtype='f8')
    dom_elev = np.full((1, nOutPolygonsHRU, ndom), hru_elev, dtype='f8')

    domType = np.full((1, nOutPolygonsHRU, ndom), 1, dtype='i4')
    if glac_dom: # NOTE, if HRU glacier area is 0, midGlce_dom should be 0
        lyrDepth[1,:,0:midToto_glac0] = lyrDepth_glac[:,:,0]
        lyrHeight[1,:,0:midToto_glac0+1] = lyrHeight_glac[:,:,0]
        midGlce_dom[:,0,1] = midGlce
        midSoil_dom[:,0,1] = midSoil_glac[:,0]
        domType[:,0,1] = 2 # glacier accumulation
        dom_area[:,0,1] = hru_area * glacAcc_frac
        dom_elev[:,0,1] = glacAcc_elev
        if clean_dom:
            lyrDepth[2,:,0:midToto_glac0] = lyrDepth_glac[:,:,0]
            lyrHeight[2,:,0:midToto_glac0+1] = lyrHeight_glac[:,:,0]
            midGlce_dom[:,0,2] = midGlce
            midSoil_dom[:,0,2] = midSoil_glac[:,0]
            domType[:,0,2] = 3 # glacier ablation clean
            dom_area[:,0,2] = hru_area * glacCln_frac
            dom_elev[:,0,2] = glacCln_elev
        if (debr_dom > 0):
            lyrDepth[indDebr,:,0:midToto_glac0] = lyrDepth_glac[:,:,1]
            lyrHeight[indDebr,:,0:midToto_glac0+1] = lyrHeight_glac[:,:,1]
            midGlce_dom[:,0,indDebr] = midGlce
            midSoil_dom[:,0,indDebr] = midSoil_glac[:,1]
            domType[:,0,indDebr] = 4 # glacier ablation debris
            dom_area[:,0,indDebr] = hru_area * glacDbr_frac
            dom_elev[:,0,indDebr] = glacDbr_elev

    if wtld_dom: # NOTE, if HRU wetland area is 0, midLake_dom should be 0
        lyrDepth[indWtld,:,0:midToto_wtld] = lyrDepth_wtld
        lyrHeight[indWtld,:,0:midToto_wtld+1] = lyrHeight_wtld
        midLake_dom[:,0,indWtld] = midLake
        midSoil_dom[:,0,indWtld] = midSoil_wtld
        domType[:,0,indWtld] = 5
        wtld_elev = hru_elev  # assume wetland elev same as upland
        dom_area[:,0,indWtld] = hru_area * wtld_frac 
        dom_elev[:,0,indWtld] = wtld_elev           # assume wetland elev same as upland

    dom_area[:,0,0] = hru_area * upld_frac
    dom_elev[:,0,0] = upld_elev

    lyrDepth = lyrDepth.transpose()
    lyrHeight = lyrHeight.transpose()

    # initialize netcdf file by storing dimensions and hru variable
    nc_out = writeNC_dims(nc_out_name, scalarv, midSoil, midToto, ifcToto,
                        outPolyIDs, gruIDs, hru_type, ndom, glac)

    # === now loop through variables and write
    #  this could be done by looping through the input state file and xferring values
    #  domain order has to be upland, glacier ablation, glacier accumulation, wetland
    #   has to have upland, then if has glacier has to have ablation and accumulation

    # domType
    writeNC_state_vars_HRU_DOM(nc_out, 'domType', 'scalarv', 'f8', domType)

    # layer Depth, Height
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerDepth', 'midToto', 'f8', lyrDepth)        # Depth
    writeNC_state_vars_HRU_DOM(nc_out, 'iLayerHeight', 'ifcToto', 'f8', lyrHeight)      # Height

    # nSoil, nSnow, nGlce, nLake
    writeNC_state_vars_HRU_DOM(nc_out, 'nSoil', 'scalarv', 'f8', midSoil_dom)           # nSoil
    writeNC_state_vars_HRU_DOM(nc_out, 'nSnow', 'scalarv', 'f8', scalar0 )              # nSnow start at 0
    writeNC_state_vars_HRU_DOM(nc_out, 'nGlce', 'scalarv', 'f8', midGlce_dom)             # nGlce
    writeNC_state_vars_HRU_DOM(nc_out, 'nLake', 'scalarv', 'f8', midLake_dom)           # nLake

    # dT
    newVarVals = np.full((ndom, 1, nOutPolygonsHRU), dT, dtype='f8')
    writeNC_state_vars_HRU_DOM(nc_out, 'dt_init', 'scalarv', 'f8', newVarVals)

    # area and elevation
    writeNC_state_vars_HRU_DOM(nc_out, 'DOMarea', 'scalarv', 'f8', dom_area)            # DOMarea
    writeNC_state_vars_HRU_DOM(nc_out, 'DOMelev', 'scalarv', 'f8', dom_elev)            # DOMelev

    # SWE, SnowDepth, SfcMeltPond, SnowAlbedo, CanopyLiq, CanopyIce, glacMass4AreaChange
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSWE', 'scalarv', 'f8', scalar0)           # SWE
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSnowDepth', 'scalarv', 'f8', scalar0)     # SnowDepth
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSfcMeltPond', 'scalarv', 'f8', scalar0)   # SfcMeltPond
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarSnowAlbedo', 'scalarv', 'f8', scalar0)    # SnowAlbedo
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyLiq', 'scalarv', 'f8', scalar0)     # CanopyLiq
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyIce', 'scalarv', 'f8', scalar0)     # CanopyIce
    writeNC_state_vars_HRU_DOM(nc_out, 'glacMass4AreaChange','scalarv', 'f8', scalar0)  # glacMass4AreaChange

    # glacier area
    writeNC_state_vars_GRU_VEC(nc_out, 'glacAblArea', 'glac', 'f8', ablArea)  
    writeNC_state_vars_GRU_VEC(nc_out, 'glacAccArea', 'glac', 'f8', accArea)
    writeNC_state_vars_GRU_VEC(nc_out, 'glacId', 'glac', 'i8', gridIDs)

    # glacier volume
    writeNC_state_vars_GRU(nc_out, 'basin__GlacierStorage', 'f8', totVolume*0.9167)     # GlacierStorage in Gt
 
    # scalar CanairTemp, CanopyTemp, AquiferStorage
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanairTemp', 'scalarv', 'f8', scalar283)  # CanairTemp, does not exist for glacier or lake
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarCanopyTemp', 'scalarv', 'f8', scalar283)  # CanopyTemp, does not exist for glacier or lake
    writeNC_state_vars_HRU_DOM(nc_out, 'scalarAquiferStorage', 'scalarv', 'f8',scalar1) # AquiferStorage, does not exist for glacier

    # layer MatricHead, Temp, VolFracLiq, VolFracIce
    if glac_dom:
        for k,h in enumerate(outPolyIDs):
            for i in range(ndom): # put ice layers at -5 C and all ice similar to Giese et al. 2020, otherwise need to spin up 40 yrs
                if i>0 and i<4: 
                    if (domType[k,0,i] == 2) or (domType[k,0,i] == 3): # glacier accumulation or ablation clean
                        ind = 0
                    elif (domType[k,0,i] == 4): # glacier ablation debris
                        ind = 1
                    toto283[   midSoil_glac[k,ind]:(midGlce+midSoil_glac[k,ind]),:,i] = 268.16 # or 273.16?
                    toto0[     midSoil_glac[k,ind] ,:,i] = 0.90 # could be 0.9 per Bradford et al. 2009, less air as go deeper
                    toto0[     midSoil_glac[k,ind]+1,:,i] = 0.91
                    toto0[     midSoil_glac[k,ind]+2,:,i] = 0.93
                    toto0[     midSoil_glac[k,ind]+3,:,i] = 0.95 
                    toto0[     midSoil_glac[k,ind]+4,:,i] = 0.98
                    totopoint2[midSoil_glac[k,ind]:(midGlce+midSoil_glac[k,ind]),:,i] = 0.0

    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerMatricHead', 'midSoil', 'f8', soilneg1)   # MatricHead
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerTemp', 'midToto', 'f8', toto283)          # Temp
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerVolFracLiq', 'midToto', 'f8', totopoint2) # VolFracLiq
    writeNC_state_vars_HRU_DOM(nc_out, 'mLayerVolFracIce', 'midToto', 'f8', toto0)      # VolFracIce


    nc_out.close()


