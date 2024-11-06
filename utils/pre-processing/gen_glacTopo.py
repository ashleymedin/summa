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
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim1,newVarDim2,'ngl','gru',),fill_value='-999')    
    else:
        ncvar = nc_out.createVariable(newVarName, newVarType, (newVarDim1,newVarDim2,'ngl','gru',),fill_value='-999.0')    
    ncvar[:,:] = newVarVals   # store data in netcdf file

# write dimensions and dimension variables to netcdf output file
def writeNC_dims(fn, grus, hru_type, ngl, Ny, Nx):    
    """ Write <vars> array in netCDF4 file,<fn> and variable of
        <varname> """
    print("writing output file")
    nc_out = nc4.Dataset(fn, 'w', format='NETCDF4')

    # Create dimensions
    dim_gru = nc_out.createDimension('gru', len(grus))
    dim_ngl = nc_out.createDimension('ngl', ngl) # max number of glaciers in any GRU
    dim_ny = nc_out.createDimension('nygrid', Ny) # max number of cells in glacier bed
    dim_nx = nc_out.createDimension('nxgrid', Nx) # max number of cells in glacier bed

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

def make_base_elevation(height_file,widehead_use,dx,dy,max_elev,min_elev,x_max,y_max):
    """ Create a base elevation array for each glacier in a GRU"""
    #x_max = 200*100. # meters, max length of the glacier bed
    #max_elev = 3000. # meters, max elevation of glacier bed (at x=0)
    #min_elev = 1000. # meters, min elevation of glacier bed (at x=x_max)
    #y_max = 400. # meters, width of glacier bed at terminus
    #DEM grid spacing, best if evenly divides x_max and y_max
    #dx = 100.        # DEM grid spacing down glacier, default 100m
    #dy = 50.         # DEM grid spacing across glacier, default 25m

    # Pick Bed Shape, choose constant, cliff, parabolic, trapezoid
    bed_shape = "parabolic"

    # Choose if how far from the head of main stem glacier the y_max width should be different than terminus, can be 0 
    # Since there are currently no boundaries, just high walls, this is unstable
    if widehead_use:
        widehead = 50.*dx # meters, 0 for no change
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
    trib_x_max = np.ones(n_trib)*5*dx     # meters, max length of the tributary
    trib_max_elev = np.ones(n_trib)*3000. # meters, max elevation of tributary
    trib_min_elev = np.ones(n_trib)*2800. # meters, min elevation of tributary (where joins the glacier)
    trib_y_max = np.ones(n_trib)*200.     # meters, width of tributary at terminus

    # glacier side tributary joins the main stem 
    # -1 is left looking toward the terminus
    # +1 is right looking toward the terminus
    trib_join = np.ones(n_trib)*1 

    # Can change these bed_shape specific values, some need to be fixed to match OGGM
    # Can change these bed_shape specific values, some need to be fixed to match OGGM
    if bed_shape=="constant": 
        y_max = y_max      # constant bed parameters set in every bed_shape, set above
    elif bed_shape=="cliff": 
        cliff_hgt = 250.   # cliff drop (meters)
        upcliff = 30.*dx   # cliff distance down glacier from head (meters)
    elif bed_shape=="parabolic": 
        shape = 5.e-3      # parabola shape, default 5.e-3 which is flattish, if use 5.e-2 very curved
    elif bed_shape=="trapezoid": 
        lambdas = 0.25     # default 2, or atan(2/lambdas)=angle of side wall from horizontal, 2 is 45 deg, 0.25 very steep
        bot_wid = y_max/2. # bottom width of trapezoid (m)
    else:
        raise ValueError("Invalid bed_shape: {}".format(bed_shape))

    # Can change these and input to OGGM, can test tributaries with widehead>0. and n_trib>=0
    if widehead>0.:
        head_y_max = y_max*2. # (meters) width at head, could be smaller than y_max
    else:
        head_y_max = y_max # set these equal if no widehead

    # Make the grid
    Nx = int(np.ceil(x_max/dx))
    dly=dx*dy
    x = np.arange(0.5, Nx + 0.5, 1)*dx
    dist = x

    Ny = int(np.ceil(y_max/dy))
    if widehead>0.: Ny = int(np.ceil(head_y_max/dy))
    y = np.arange(0.5, Ny+0.5, 1)*dy

    # Bed array initialize
    B = np.zeros((Ny,Nx)) 
    glacierMask = np.ones((Ny, Nx), dtype=bool) # glacier mask of True values
    xm, ym = np.meshgrid(x, y) # meshgrid of x and y

    # elevation[m] of constant slope bedrock at the flowline grid points
    slope_h = np.linspace(max_elev,min_elev,Nx) 

    # specify the bed elevation at each grid point
    if bed_shape == "constant" or bed_shape == "cliff":
        B = np.tile(slope_h, (Ny, 1))
        if bed_shape == "cliff":
            icliff = int(upcliff/dx)
            B[:, icliff:] -= cliff_hgt

    elif bed_shape == "trapezoid":
        y0_max = np.maximum(head_y_max, y_max)
        base = np.zeros_like(y)
        base[y <= (y0_max-bot_wid)/2.] = (y0_max-bot_wid)/2. - y[y <= (y0_max-bot_wid)/2.]
        base[y >= (y0_max-bot_wid)/2. + bot_wid] = -( (y0_max+bot_wid)/2.-y[y >= (y0_max-bot_wid)/2. + bot_wid] )
        B = base[:, np.newaxis]*(2./lambdas) + slope_h[np.newaxis, :]

    elif bed_shape == "parabolic":
        y0_max = np.maximum(head_y_max, y_max)
        base = np.abs(y - y0_max/2.)
        B = shape*(base[:,np.newaxis]**2.) + slope_h[np.newaxis,:]
        
    y_flow = int(Ny/2.) #if odd Ny, this will be center FD coordinate, otherwise will not be center
    #bed_flow = B[y_flow]

    # cut off bed if out of bounds for plotting and for MB calculations
    y0_max = np.maximum(head_y_max, y_max)
    center = y[y_flow] # center of glacier grid, needs to be in the middle of a grid cell

    # cut off the bed outside of width
    wid_max_xm = head_y_max * (xm <= widehead) + y_max * (xm > widehead)
    glacierMask[(ym > center + wid_max_xm/2.) | (ym < center - wid_max_xm/2.)] = False

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

    return Ny, Nx, B, S, glacierMask, area, volume


############################################
#                  Main                    #
############################################

if __name__ == '__main__':
    if testing:
        # hardwired for testing
        nc_attribute_name = 'attributes.nc'
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
        hru_type = args[1]              # 'int' or 'string'
 
    # hardwired to forcing formats (hru index rather than grid)
    outPolyIDs, hru_elev, hru_area, gruIDs, hru2gru = getOutputPolyIDs(nc_attribute_name)        
    nOutPolygonsHRU = len(outPolyIDs)
    nOutPolygonsGRU = len(gruIDs)

    nc_out_name = nc_attribute_name[:-3] + '_base_topo_glac.nc'
    print("created glacier base topo file")

    nc_outS_name = nc_attribute_name[:-3] + '_surf_topo_glac.nc'
    print("created glacier base topo file")


    # hardwired for testing, get from RGI2000-v7.0-G overlay on GRUs
    ngl = 2 # max glaciers in any GRU
    dim_ny = 10 # max number of cells in glacier bed across glacier
    dim_nx = 200 # max number of cells in glacier bed down glacier
    # initialize netcdf file by storing dimensions and hru variable
    widehead_use = True
    if widehead_use: 
        nc_out = writeNC_dims(nc_out_name, gruIDs, hru_type, ngl, 2*dim_ny, dim_nx)
        nc_outS = writeNC_dims(nc_outS_name, gruIDs, hru_type, ngl, 2*dim_ny, dim_nx)
    else:
        nc_out = writeNC_dims(nc_out_name, gruIDs, hru_type, ngl, dim_ny, dim_nx)
        nc_outS = writeNC_dims(nc_outS_name, gruIDs, hru_type, ngl, dim_ny, dim_nx)    
    height_file = '/Users/amedin/Research/USask/GlacierPython/my_code/mymodel_height_linear100.npz' # 'none' if no surface height file
    totArea = np.zeros(nOutPolygonsGRU, dtype='f8')
    nGlacier = np.zeros(nOutPolygonsGRU, dtype='i4')
    glacArea = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    glacVol = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    ela = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    Ny = np.zeros((1, ngl, nOutPolygonsGRU), dtype='i4')
    Nx = np.zeros((1, ngl, nOutPolygonsGRU), dtype='i4')
    dy = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    dx = np.zeros((1, ngl, nOutPolygonsGRU), dtype='f8')
    if widehead_use: 
        B = np.zeros((2*dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
        S = np.zeros((2*dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
        glacierMask = np.zeros((2*dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='i4')
        cell2hru = np.zeros((2*dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
    else:
        B = np.zeros((dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
        S = np.zeros((dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
        glacierMask = np.zeros((dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='i4')
        cell2hru = np.zeros((dim_ny, dim_nx, ngl, nOutPolygonsGRU), dtype='f8')
    totElev = np.zeros(nOutPolygonsGRU)
    np.random.seed(13) # for reproducibility
    for i,g in enumerate(gruIDs):
        totArea[i] = hru_area[hru2gru==g].sum()
        totElev[i] = hru_elev[hru2gru==g].mean() # really only makes sense if all of glacier is in one HRU
        nGlacier[i] =2 # number of glaciers in each GRU, should get from RGI2000-v7.0-G
        glac_frac = 0.5 # fraction of glacier area in each GRU, for testing
        fac = np.sqrt(totArea[i]*glac_frac/((129.1+47.7)*1e6))
        for j in range(nGlacier[i]):
            # Should get the glacier area from RGI2000-v7.0-G area_km2
            if j==0:
                # Should get these from RGI2000-v7.0-G zmax_m, zmin_m, zmed_m, or from RGI2000-v7.0-G-01_alaska-attributes.csv
                #   these numbers are from Wernike Glacier, Alaska
                area = 129.1 # km^2, glacier area
                max_elev = 2776.3 # meters, max elevation of glacier (bed at x=0) 
                min_elev = 136.8 # meters, min elevation of glacier (bed at x=x_max)
                med_elev = 1302.1 # meters, median elevation of glacier bed
                # should get the glacier thickness from Millan et al. 2022, and surface from DEM
                # then, calculate bed elevation and a grid over it
                # for now, just make a simple glacier bed
                x_max = 35366.*fac # meters, max length of the glacier bed
                y_max = 3000.*fac # meters, width of glacier bed at terminus
                #DEM grid spacing, best if evenly divides x_max and y_max
                dx0 = x_max/dim_nx        # DEM grid spacing down glacier, default 100m
                dy0 = y_max/dim_ny        # DEM grid spacing across glacier, default 25m
            if j>=1:
                #   these numbers are from Van Cleve Glacier, Alaska
                area = 47.7 # km^2, glacier area
                max_elev = 2290.1 # meters, max elevation of glacier (bed at x=0)
                min_elev =309.7 # meters, min elevation of glacier (bed at x=x_max)
                med_elev = 1232.0 # meters, median elevation of glacier bed
                x_max = 19715.*fac # meters, max length of the glacier bed
                y_max = 2000.*fac # meters, width of glacier bed at terminus
                dx0 = x_max/dim_nx        # DEM grid spacing down glacier, default 100m
                dy0 = y_max/dim_ny        # DEM grid spacing across glacier, default 25m
            
            Ny0, Nx0, B0, S0, glacierMask0, area, volume = make_base_elevation(height_file,widehead_use,dx0,dy0,max_elev,min_elev,x_max,y_max)

            # need to define each grid cell to an hruID "cell2hru", 
            # grid cells have vars X,Y coordinates, B, S, glacierMask, and cell2hru, (dim glacierID, dim gruID) xxx[:nCells,j,i]
            ela[0,j,i] = med_elev # meters, surface elevation at ELA
            glacArea[0,j,i] = area
            glacVol[0,j,i] = volume
            Ny[0,j,i] = Ny0
            Nx[0,j,i] = Nx0
            dx[0,j,i] = dx0
            dy[0,j,i] = dy0
            B[:Ny0,:Nx0,j,i] = B0
            S[:Ny0,:Nx0,j,i] = S0
            glacierMask[:Ny0,:Nx0,j,i] = np.where(glacierMask0, 1, glacierMask[:Ny0, :Nx0, j, i])
            if nOutPolygonsGRU==nOutPolygonsHRU:
                cell2hru[:Ny0,:Nx0,j,i] = hru2gru[hru2gru==g] 
            else:
                sys.exit("ERROR, nOutPolygonsGRU must equal nOutPolygonsHRU or need to define cell2hru")

    writeNC_state_vars_GRU(nc_out, 'nGlacier', 'i4', nGlacier)
    writeNC_state_vars_GRU_VEC(nc_out, 'Ny', 'ngl', 'i4', Ny)
    writeNC_state_vars_GRU_VEC(nc_out, 'Nx', 'ngl', 'i4', Nx)
    writeNC_state_vars_GRU_VEC(nc_out, 'dx', 'ngl', 'f8', dx)
    writeNC_state_vars_GRU_VEC(nc_out, 'dy', 'ngl', 'f8', dy)
    writeNC_state_vars_GRU_GRID(nc_out, 'bed_elev_m', 'nygrid','nxgrid','f8', B)
    writeNC_state_vars_GRU_GRID(nc_out, 'glacierMask', 'nygrid','nxgrid','i4', glacierMask)
    writeNC_state_vars_GRU_GRID(nc_out, 'cell2hru', 'nygrid','nxgrid','f8', cell2hru)

    # cells with S above median elevation are accumulation area, below are ablation area
    writeNC_state_vars_GRU_VEC(nc_outS, 'glacArea_km2', 'ngl', 'f8', glacArea)
    writeNC_state_vars_GRU_VEC(nc_outS, 'glacVol_km3', 'ngl', 'f8', glacVol)
    writeNC_state_vars_GRU_VEC(nc_outS, 'ELA_elev_m', 'ngl', 'f8', ela)
    writeNC_state_vars_GRU_GRID(nc_outS, 'surface_elev_m', 'nygrid','nxgrid','f8', S)

    nc_out.close()
    nc_outS.close()