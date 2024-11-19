#!/usr/bin/env python
''' Create a glacier topographic file for SUMMA'''
#
# Author:  Ashley Van Beusekom, 2024
#
# Notes: make a simple glacier topographic file with input glacier area
#  median,min,max elevation, at least one glacier thickness, and assume parabolic shape
##
# Requirements:  run with a python (eg miniconda) 2.7 that includes netCDF4
# Run as: python gen_glacTopo.py <existing_attributeFile_with_hruId> <hru_type (int|int64|str)>
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

# write gru nonscalar variables to netcdf output file
def writeNC_state_vars_GRU_VEC(nc_out, newVarName, newVarDim, newVarType, newVarVals):

    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding GRU_VEC data")
    if newVarType=='i4' or newVarType=='i8':
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'gru',),fill_value='-999')  
    else:
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim,'gru',),fill_value='-999.0')  
    ncvar[:] = newVarVals   # store data in netcdf file

# write gru grid variables to netcdf output file
def writeNC_state_vars_GRU_GRID(nc_out, newVarName, newVarDim1,newVarDim2, newVarType, newVarVals):

    """ Write <vars>[gru] array in netCDF4 file,<fn> and variable of
        <varname> """

    print("adding GRU_GRID data")
    if newVarType=='i4' or newVarType=='i8':
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim1,newVarDim2,'grid','gru',),fill_value='-999')    
    else:
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim1,newVarDim2,'grid','gru',),fill_value='-999.0')    
    ncvar[:,:] = newVarVals   # store data in netcdf file

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn, grus, hru_type, glac, nx, ny):    
    """ Write <vars> array in netCDF4 file,<fn> and variable of
        <varname> """
    print("writing output file")
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_gru = nc_out.createDimension('gru', len(grus))
    dim_ngl = nc_out.createDimension('grid', glac) # max number of glaciers in any GRU
    dim_nx = nc_out.createDimension('xgrid', nx) # max number of cells in glacier bed
    dim_ny = nc_out.createDimension('ygrid', ny) # max number of cells in glacier bed

    # --- Create HRU ID variable (can be either int or string)
    if hru_type == 'str':
        # string HRU (need to add string length)
        max_strlen = 20  # EC
        dim_str = nc_out.createDimension('strlen', max_strlen)
        gruId = nc_out.createVariable('gruId', 'S1', ('gru', 'strlen'),fill_value='-999')
        gruId[:] = nc4.stringtochar(np.asarray(np.unique(grus),
                                  dtype='S{}'.format(max_strlen)))

    elif hru_type == 'int64':
        # integer HRU
        gruId = nc_out.createVariable('gruId', 'i8', ('gru', ),fill_value='-999')
        gruId[:] = grus

    elif hru_type == 'int':
        # integer HRU
        gruId = nc_out.createVariable('gruId', 'int', ('gru', ),fill_value='-999')
        gruId[:] = grus

    else:
        # not recognized
        sys.exit("ERROR, hru_type not recognized: must be str, int64, or int")

    # add attribute    
    gruId.long_name = 'GRU ID'

    return nc_out
    # leave netcdf file open

def make_base_elevation(height_file,widehead_use,dy,dx,max_elev,min_elev,y_max,x_max):
    """ Create a base elevation array for each glacier in a GRU"""
    #y_max = 200*100. # meters, max length of the glacier bed
    #max_elev = 3000. # meters, max elevation of glacier bed (at y=0)
    #min_elev = 1000. # meters, min elevation of glacier bed (at y=y_max)
    #x_max = 400. # meters, width of glacier bed at terminus
    #DEM grid spacing, best if evenly divides y_max and x_max
    #dy = 100.        # DEM grid spacing down glacier, default 100m
    #dx = 50.         # DEM grid spacing across glacier, default 25m

    # Pick Bed Shape, choose constant, cliff, parabolic, trapezoid
    bed_shape = "parabolic"

    # Choose if how far from the head of main stem glacier the x_max width should be different than terminus, can be 0 
    # Since there are currently no boundaries, just high walls, this is unstable
    if widehead_use:
        widehead = 50.*dy # meters, 0 for no change
    else:
        widehead = 0.

    # Choose how many tributaries to add, can be 0
    # Tributaries position and order can be changed in the code
    # Note, for realism will need to define flowline structure including order
    # It makes sense to have the mass balance the same throughout all parts of the glacier
    #   and the bed shape, but the shape parameters should be allowed to change
    # Since there are currently no boundaries, just high walls, this is unstable
    n_trib = 0

    # give a value per tributary, will not be used if n_trib==0
    trib_y_max = np.ones(n_trib)*5*dy     # meters, max length of the tributary
    trib_max_elev = np.ones(n_trib)*3000. # meters, max elevation of tributary
    trib_min_elev = np.ones(n_trib)*2800. # meters, min elevation of tributary (where joins the glacier)
    trib_x_max = np.ones(n_trib)*200.     # meters, width of tributary at terminus

    # glacier side tributary joins the main stem 
    # -1 is left looking toward the terminus
    # +1 is right looking toward the terminus
    trib_join = np.ones(n_trib)*1 

    # Can change these bed_shape specific values, some need to be fixed to match OGGM
    # Can change these bed_shape specific values, some need to be fixed to match OGGM
    if bed_shape=="constant": 
        x_max = x_max      # constant bed parameters set in every bed_shape, set above
    elif bed_shape=="cliff": 
        cliff_hgt = 250.   # cliff drop (meters)
        upcliff = 30.*dy   # cliff distance down glacier from head (meters)
    elif bed_shape=="parabolic": 
        shape = 5.e-3      # parabola shape, default 5.e-3 which is flattish, if use 5.e-2 very curved
    elif bed_shape=="trapezoid": 
        lambdas = 0.25     # default 2, or atan(2/lambdas)=angle of side wall from horizontal, 2 is 45 deg, 0.25 very steep
        bot_wid = x_max/2. # bottom width of trapezoid (m)
    else:
        raise ValueError("Invalid bed_shape: {}".format(bed_shape))

    # Can change these and input to OGGM, can test tributaries with widehead>0. and n_trib>=0
    if widehead>0.:
        head_x_max = x_max*2. # (meters) width at head, could be smaller than x_max
    else:
        head_x_max = x_max # set these equal if no widehead

    # Make the grid
    ny = int(np.ceil(y_max/dy))
    dly=dy*dx
    y = np.arange(0.5, ny + 0.5, 1)*dy

    nx = int(np.ceil(x_max/dx))
    if widehead>0.: nx = int(np.ceil(head_x_max/dx))
    x = np.arange(0.5, nx+0.5, 1)*dx

    # Bed array initialize
    B = np.zeros((nx,ny)) 
    glacierMask = np.ones((nx, ny), dtype=bool) # glacier mask of True values
    ym, xm = np.meshgrid(y, x) # meshgrid of y and x

    # elevation[m] of constant slope bedrock at the flowline grid points
    slope_h = np.linspace(max_elev,min_elev,ny) 

    # specify the bed elevation at each grid point
    if bed_shape == "constant" or bed_shape == "cliff":
        B = np.tile(slope_h, (nx, 1))
        if bed_shape == "cliff":
            icliff = int(upcliff/dy)
            B[:, icliff:] -= cliff_hgt

    elif bed_shape == "trapezoid":
        x0_max = np.maximum(head_x_max, x_max)
        base = np.zeros_like(x)
        base[x <= (x0_max-bot_wid)/2.] = (x0_max-bot_wid)/2. - x[x <= (x0_max-bot_wid)/2.]
        base[x >= (x0_max-bot_wid)/2. + bot_wid] = -( (x0_max+bot_wid)/2.-x[x >= (x0_max-bot_wid)/2. + bot_wid] )
        B = base[:, np.newaxis]*(2./lambdas) + slope_h[np.newaxis, :]

    elif bed_shape == "parabolic":
        x0_max = np.maximum(head_x_max, x_max)
        base = np.abs(x - x0_max/2.)
        B = shape*(base[:,np.newaxis]**2.) + slope_h[np.newaxis,:]
        
    x_flow = int(nx/2.) #if odd nx, this will be center FD coordinate, otherwise will not be center
    #bed_flow = B[x_flow]

    # cut off bed if out of bounds for plotting and for MB calculations
    x0_max = np.maximum(head_x_max, x_max)
    center = x[x_flow] # center of glacier grid, needs to be in the middle of a grid cell

    # cut off the bed outside of width
    wid_max_ym = head_x_max * (ym <= widehead) + x_max * (ym > widehead)
    glacierMask[(xm > center + wid_max_ym/2.) | (xm < center - wid_max_ym/2.)] = False

    # initialize to make sure it doesn't grow out of bounds
    B = np.where(glacierMask, B, B + 1000)

    # surface elevation
    if height_file=='none':
        S = np.copy(B)+ np.where(glacierMask, 1, 0) # put 1 m of ice on top of bed
    else:
        npzfile=np.load(height_file)
        H = npzfile['height']
        S = B + H

    H = S - B # ice thickness
    area = np.sum(np.where(H>0, dly, 0))*1.e-6 # glacier area km2
    volume = np.sum(np.where(H>0, H*dly, 0))*1.e-9 # glacier volume km3

    return nx, ny, B, S, glacierMask, area, volume


############################################
#                  Main                    #
############################################

if __name__ == '__main__':
    if testing:
        # hardwired for testing
        nc_attribute_name = 'attributes.nc'
        hru_type = 'int'
        nc_out_name0 = 'coldState.nc'

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
        nc_out_name0 = args[1]          # output cold-state file
        hru_type = args[2]              # 'int' or 'string'
 
    # hardwired to forcing formats (hru index rather than grid)
    outPolyIDs, hru_elev, hru_area, gruIDs, hru2gru = getOutputPolyIDs(nc_attribute_name)        
    nOutPolygonsHRU = len(outPolyIDs)
    nOutPolygonsGRU = len(gruIDs)

    nc_out_name = nc_attribute_name[:-3] + '_glacBedTopo.nc'
    print("created glacier bed topo file")

    nc_outS_name = nc_out_name0[:-3] + '_glacSurfTopo.nc'
    print("created glacier surface topo file")

    # hardwired for testing, get from RGI2000-v7.0-G overlay on GRUs
    glac = 2 # max glaciers in any GRU
    dim_nx = 10 # max number of cells in glacier bed across glacier
    dim_ny = 200 # max number of cells in glacier bed down glacier
    # initialize netcdf file by storing dimensions and hru variable
    widehead_use = True
    if widehead_use: 
        nc_out = writeNC_dims(nc_out_name, gruIDs, hru_type, glac, 2*dim_nx, dim_ny)
        nc_outS = writeNC_dims(nc_outS_name, gruIDs, hru_type, glac, 2*dim_nx, dim_ny)
    else:
        nc_out = writeNC_dims(nc_out_name, gruIDs, hru_type, glac, dim_nx, dim_ny)
        nc_outS = writeNC_dims(nc_outS_name, gruIDs, hru_type, glac, dim_nx, dim_ny)    
    height_file = '/Users/amedin/Research/USask/GlacierPython/my_code/mymodel_height_linear100.npz' # 'none' if no surface height file
    totArea = np.zeros(nOutPolygonsGRU, dtype='f8')
    nGlacier = np.zeros(nOutPolygonsGRU, dtype='i4')
    gridId = np.zeros((1, glac, nOutPolygonsGRU), dtype='i8')
    nx = np.zeros((1, glac, nOutPolygonsGRU), dtype='i4')
    ny = np.zeros((1, glac, nOutPolygonsGRU), dtype='i4')
    dx = np.zeros((1, glac, nOutPolygonsGRU), dtype='f8')
    dy = np.zeros((1, glac, nOutPolygonsGRU), dtype='f8')
    debris_thick = np.zeros((1, glac, nOutPolygonsGRU), dtype='f8')
    AAR_frac = np.zeros((1, glac, nOutPolygonsGRU), dtype='f8')
    stage_frac = np.zeros((1, glac, nOutPolygonsGRU), dtype='f8')

    if widehead_use: 
        B = np.zeros((dim_ny, 2*dim_nx, glac, nOutPolygonsGRU), dtype='f8')
        S = np.zeros((dim_ny, 2*dim_nx, glac, nOutPolygonsGRU), dtype='f8')
        glacierMask = np.zeros((dim_ny, 2*dim_nx, glac, nOutPolygonsGRU), dtype='i4')
        cell2hruId  = np.zeros((dim_ny, 2*dim_nx, glac, nOutPolygonsGRU), dtype='i8')
    else:
        B = np.zeros((dim_ny, dim_nx, glac, nOutPolygonsGRU), dtype='f8')
        S = np.zeros((dim_ny, dim_nx, glac, nOutPolygonsGRU), dtype='f8')
        glacierMask = np.zeros((dim_ny, dim_nx, glac, nOutPolygonsGRU), dtype='i4')
        cell2hruId  = np.zeros((dim_ny, dim_nx, glac, nOutPolygonsGRU), dtype='i8')
    totElev = np.zeros(nOutPolygonsGRU)
    np.random.seed(13) # for reproducibility
    for i,g in enumerate(gruIDs):
        totArea[i] = hru_area[hru2gru==g].sum()
        totElev[i] = hru_elev[hru2gru==g].mean() # really only makes sense if all of glacier is in one HRU
        nGlacier[i] =2 # number of glaciers in each GRU, should get from RGI2000-v7.0-G
        glac_frac = 0.5 # fraction of glacier area in each GRU, for testing
        fac = np.sqrt(totArea[i]*glac_frac/((129.1+47.7)*1e6))

        # Should get these from RGI2000-v7.0-G RGI_ID, zmax_m, zmin_m from RGI2000-v7.0-G-01_alaska-attributes.csv
        # Can get S from DEM and B from tif files here https://www.sedoo.fr/theia-publication-products/?uuid=55acbdd5-3982-4eac-89b2-46703557938c
        # Should have two ablation zones, one with debri (~30% often in alaska) and one without
        # Can get AAR and stage(deb_area/abl_area) or deb_km2, abl_km2, gl_km2 by RGI ID in Pelicotti product, if >2km2 glacier, and just debri maps if 1-2 km2 ... would have to calculate areas (only calculated by region) https://zenodo.org/records/3866466
	    #   - How would you recalculate debri area and mean elevation as glacier melted? keep it always the same percentage?? (==stage)
	    #   - What is thickness? Rounce has thickness tif files here https://nsidc.org/data/hma_dte/versions/1#anchor-data-access-tools
        #   - maybe leave thickness and percent of HRU domain ablation area as a constant in attributes per glacier, per HRU?

        for j in range(nGlacier[i]):
            if j==0:
                #   these numbers are from Wernike Glacier, Alaska
                RGI_ID = 'RGI60-01.00004'
                gridId0 = int(RGI_ID.split('-')[1].split('.')[0]+RGI_ID.split('-')[1].split('.')[1])
                #area_km2 = 129.1 # km^2, glacier area
                max_elev = 2776.3 # meters, max elevation of glacier (bed at y=0) 
                min_elev = 136.8 # meters, min elevation of glacier (bed at y=y_max)
                # should get the glacier thickness from Millan et al. 2022, and surface from DEM
                # then, calculate bed elevation and a grid over it
                # for now, just make a simple glacier bed
                y_max = 35366.*fac # meters, max length of the glacier bed
                x_max = 3000.*fac # meters, width of glacier bed at terminus
                #DEM grid spacing, best if evenly divides y_max and x_max
                dy0 = y_max/dim_ny        # DEM grid spacing down glacier, default 100m
                dx0 = x_max/dim_nx        # DEM grid spacing across glacier, default 25m
                debris_m = 0.25 # if debris covered, need to add a debris thickness
                AAR = 0.5 # acculation area ratio, percent of glacier area that is accumulation area
                stage = 0.3 # percent of ablation area that is debris covered
            if j>=1:
                #   these numbers are from Van Cleve Glacier, Alaska
                RGI_ID = 'RGI60-01.00005'
                gridId0 = int(RGI_ID.split('-')[1].split('.')[0]+RGI_ID.split('-')[1].split('.')[1])
                #area_km2 = 47.7 # km^2, glacier area
                max_elev = 2290.1 # meters, max elevation of glacier (bed at y=0)
                min_elev =309.7 # meters, min elevation of glacier (bed at y=y_max)
                y_max = 19715.*fac # meters, max length of the glacier bed
                x_max = 2000.*fac # meters, width of glacier bed at terminus
                dy0 = y_max/dim_ny        # DEM grid spacing down glacier, default 100m
                dx0 = x_max/dim_nx        # DEM grid spacing across glacier, default 25m
                debris_m = 0.15 # if debris covered, need to add a debris thickness
                AAR = 0.5 # acculation area ratiox, percent of glacier area that is accumulation area
                stage = 0.3 # percent of ablation area that is debris covered
            
            nx0, ny0, B0, S0, glacierMask0, area_km2, volume = make_base_elevation(height_file,widehead_use,dy0,dx0,max_elev,min_elev,y_max,x_max)

            # need to define each grid cell to an hruID "cell2hruId", 
            # grid cells have vars X,Y coordinates, B, S, glacierMask, and cell2hruId, (dim glacierID, dim gruID) xxx[:nCells,j,i]
            gridId[0,j,i] = gridId0
            nx[0,j,i] = nx0
            ny[0,j,i] = ny0
            dy[0,j,i] = dy0
            dx[0,j,i] = dx0
            B[:ny0,:nx0,j,i] = B0.transpose()
            S[:ny0,:nx0,j,i] = S0.transpose()
            debris_thick[0,j,i] = debris_m
            AAR_frac[0,j,i] = AAR
            stage_frac[0,j,i] = stage
            glacierMask[:ny0,:nx0,j,i] = np.where(glacierMask0.transpose(), 1, glacierMask[:ny0,:nx0, j, i])
            if nOutPolygonsGRU==nOutPolygonsHRU:
                cell2hruId[:ny0,:nx0,j,i] = hru2gru[hru2gru==g] 
            else:
                sys.exit("ERROR, nOutPolygonsGRU must equal nOutPolygonsHRU or need to define cell2hruId")

    writeNC_state_vars_GRU(nc_out, 'nGlacier', 'i4', nGlacier)
    writeNC_state_vars_GRU_VEC(nc_out, 'gridId', 'grid', 'i8', gridId)
    writeNC_state_vars_GRU_VEC(nc_out, 'nx', 'grid', 'i4', nx)
    writeNC_state_vars_GRU_VEC(nc_out, 'ny', 'grid', 'i4', ny)
    writeNC_state_vars_GRU_VEC(nc_out, 'dy', 'grid', 'f8', dy)
    writeNC_state_vars_GRU_VEC(nc_out, 'dx', 'grid', 'f8', dx)
    writeNC_state_vars_GRU_VEC(nc_out, 'debris_thick', 'grid', 'f8', debris_thick)    # not needed for SUMMA attributes, but need for layer setup
    writeNC_state_vars_GRU_VEC(nc_out, 'stage', 'grid', 'f8', stage_frac)               # not needed for SUMMA attributes, but need for domain setup
    writeNC_state_vars_GRU_GRID(nc_out, 'bed_elev','ygrid','xgrid','f8', B)
    writeNC_state_vars_GRU_GRID(nc_out, 'glacierMask','ygrid','xgrid','i4', glacierMask)
    writeNC_state_vars_GRU_GRID(nc_out, 'cell2hruId','ygrid','xgrid','i8', cell2hruId)

    # cells with S above median elevation are accumulation area, below are ablation area
    writeNC_state_vars_GRU_VEC(nc_outS, 'gridId', 'grid', 'i8', gridId)
    writeNC_state_vars_GRU_VEC(nc_outS, 'AAR', 'grid', 'f8', AAR_frac)      # not needed for SUMMA initial conditions, but need to calculate initial HRU area
    writeNC_state_vars_GRU_GRID(nc_outS, 'surface_elev','ygrid','xgrid','f8', S)

    nc_out.close()
    nc_outS.close()