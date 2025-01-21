# written by A. Van Beusekom (2023)

## Visualize statistics per GRU
## Needs:
#    SUMMA output statistics

## Special note
# SUMMA simulations have been preprocessed into single value statistics per model element, using auxiliary scripts in ~/utils
# Run:
# python hist_per_GRU.py [stat]
# where stat is rmse or maxe or kgem or rmnz or avge

# modules
import os
import matplotlib
import numpy as np
import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import copy
import pandas as pd

do_rel = False # true is plot relative to the benchmark simulation
do_hist = False # true is plot histogram instead of CDF
run_local = True # true is run on local machine, false is run on cluster
fix_units_soil = True # true is convert to storage units, only works for Soil
no_snow = False # true is only plot snow free simulations

if run_local: 
    stat = 'avge'
    viz_dir = Path('/Users/amedin/Research/USask/test_py/statistics_en')
else:
    import sys
    stat = sys.argv[1]
    viz_dir = Path('/home/avanb/scratch/statistics')
    

#method_name=['be1','sundials_1en4','be4','be8','be16','be32','sundials_1en6'] #maybe make this an argument
#plt_name=['BE1','IDAe-4','BE4','BE8','BE16','BE32','IDAe-6'] #maybe make this an argument
#method_name=['be1','be16','be32','sundials_1en6'] #maybe make this an argument
#plt_name=['BE1','BE16','BE32','SUNDIALS'] #maybe make this an argument
method_name=['be8','be8cm','be8en','sundials_1en5cm','sundials_1en5en'] 
plt_name=['BE8 common','BE8 temp','BE8 mixed','SUNDIALS temp', 'SUNDIALS enth']
#method_name=['old_be1','old_be1cm','old_be1en','be8','be8cm','be8en','sundials_1en5cm','sundials_1en5en'] 
#plt_name=['BE1 common','BE1 temp','BE1 mixed','BE8 common','BE8 temp','BE8 mixed','SUNDIALS temp', 'SUNDIALS enth']
method_name2=method_name +['sundials_1en8en']
plt_name2=plt_name +['reference soln']

num_bins = 1000

if stat == 'kgem': do_rel = False # don't plot relative to the benchmark simulation for KGE

# Define the power transformation function
def power_transform(x):
    return x ** 0.5  # Adjust the exponent as needed

# Simulation statistics file locations
use_vars = []
rep = [] # mark the repeats
use_vars = [4,4,1,1]
rep = [1,2,1,2] # mark the repeats
settings0= ['scalarSWE','scalarTotalSoilWat','scalarTotalET','scalarCanopyWat','scalarRootZoneTemp']
settings = [settings0[i] for i in use_vars]

use_vars2 = [8]
rep2 = [0] # mark the repeats
use_vars2 = [3,3]
rep2 = [1,2] # mark the repeats
settings20= ['balanceCasNrg','balanceVegNrg','balanceSnowNrg','balanceSoilNrg','balanceVegMass','balanceSnowMass','balanceSoilMass','balanceAqMass','wallClockTime']
settings2 = [settings20[i] for i in use_vars2]

viz_fil = method_name.copy()
viz_fl2 = method_name2.copy()
for i, m in enumerate(method_name):
    viz_fil[i] = m + '_hrly_diff_stats_{}.nc'
    viz_fil[i] = viz_fil[i].format(','.join(settings0))
for i, m in enumerate(method_name2):
    viz_fl2[i] = m + '_hrly_diff_bals_{}.nc'
    viz_fl2[i] = viz_fl2[i].format(','.join(['balance']))

# Specify variables of interest
plot_vars = settings.copy()
plt_titl = ['snow water equivalent','total soil water content','total evapotranspiration', 'total water on the vegetation canopy','top 3m soil temperature']
leg_titl = ['$kg~m^{-2}$', '$kg~m^{-2}$','mm~y^{-1}$','$kg~m^{-2}$','$K$']
if (len(use_vars)+len(use_vars2)>1): 
    plt_titl = [f"({chr(97+n)}) {plt_titl[i]}" for n,i in enumerate(use_vars)]
else:
    plt_titl = [f"{plt_titl[i]}" for n,i in enumerate(use_vars)]
leg_titl = [leg_titl[i] for i in use_vars]

plot_vars2 = settings2.copy()
plt_titl2 = ['canopy air space enthalpy balance','vegetation enthalpy balance','snow enthalpy balance','soil enthalpy balance','vegetation mass balance','snow mass balance','soil mass balance','aquifer mass balance', 'wall clock time']
leg_titl2 = ['$W~m^{-3}$'] * 4 + ['$kg~m^{-3}~s^{-1}$'] * 3 + ['$kg~m^{-2}~s^{-1}$']+ ['$s$']
if fix_units_soil: leg_titl2 = ['$kJ~m^{-2}$'] * 4 + ['$kg~m^{-2}'] * 4 + ['$s$']
if (len(use_vars)+len(use_vars2)>1): 
    plt_titl2 = [f"({chr(97+n + len(use_vars))}) {plt_titl2[i]}" for n,i in enumerate(use_vars2)]
else:
    plt_titl2 = [f"{plt_titl2[i]}" for n,i in enumerate(use_vars2)]
leg_titl2 = [leg_titl2[i] for i in use_vars2]

if do_hist:
    fig_fil = 'Hrly_diff_hist_{}_{}_zoom'
else:
    fig_fil = 'Hrly_diff_cdf_{}_{}_zoom'
if do_rel: fig_fil = fig_fil+'_rel'
if no_snow: fig_fil = fig_fil + '_nosnow'
fig_fil = fig_fil +'_compressed.png'
fig_fil = fig_fil.format(','.join(settings),stat)

if stat == 'avge':
    stat2 = 'mean'
    maxes = [99,7,99,99,0.28]
    if do_rel: maxes = [0.4,0.007,0.6,0.15,0.0015]
if stat == 'rmse' or stat=='rmnz':
    stat2 = 'mean'
    maxes = [2,15,250,0.08,200]
    if do_rel: maxes = [0.4,0.007,0.6,0.15,0.0015]
if stat == 'maxe':
    stat2 = 'amax'
    maxes_m = [99,15,99,99,7.5]
    if do_rel: maxes_m = [0.4,0.007,0.6,0.15,0.0015]
    if stat == 'maxe': maxes = maxes_m
if stat == 'kgem':
    stat2 = 'mean'
    maxes = [0.9,0.9,0.9,0.9,0.9]
maxes = [maxes[i] for i in use_vars]
for i in range(len(maxes)):
    #if rep[i]==2: maxes[i] = maxes[i]*2.5 #clunky way to increase the range for the second repeat
    if rep[i]==2: maxes[i] = maxes_m[use_vars[i]] #clunky way to increase the range for the second repeat

if stat2 == 'mean':
    maxes2 = [1e2,1e2,1e2,1e2]+[1e-7,1e-5,1e-7,1e-8] + [5e-2]
if stat2 == 'amax':
    maxes2 = [1e4,1e4,1e4,1e4]+[1e-5,1e-3,1e-5,1e-6] + [2.0]
maxes2 = [maxes2[i] for i in use_vars2]
for i in range(len(maxes2)):
    if rep2[i]==2: maxes2[i] = maxes2[i]*1e2 #clunky way to increase the range for the second repeat

summa = {}
summa1 = {}
if len(use_vars)>0:
    for i, m in enumerate(method_name):
        # Get the aggregated statistics of SUMMA simulations
        summa[m] = xr.open_dataset(viz_dir/viz_fil[i])
if len(use_vars2)>0:
    for i, m in enumerate(method_name2):
        summa1[m] = xr.open_dataset(viz_dir/viz_fl2[i])

if no_snow:
    summa[method_name[0]] = xr.open_dataset(viz_dir/viz_fil[0]) # will be a problem if this does not exist
    if len(use_vars)>0:
        for m in method_name:
            summa[m] = summa[m].where(summa[method_name[0]]['scalarSWE'].sel(stat='mean_ben') == 0)
    if len(use_vars2)>0:
        for m in method_name2:
            summa1[m] = summa1[m].where(summa[method_name[0]]['scalarSWE'].sel(stat='mean_ben') == 0)
    
##Figure

plt.rcParams['xtick.color'] = 'black'
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.color'] = 'black'
plt.rcParams['ytick.major.width'] = 2
# fix size for now
ncol = 4
nrow = 2

if 'compressed' in fig_fil:
    plt.rcParams.update({'font.size': 27})
else:
    plt.rcParams.update({'font.size': 100})

if 'compressed' in fig_fil:
    fig,axs = plt.subplots(nrow,ncol,figsize=(17*ncol,17*nrow))
else:
    fig,axs = plt.subplots(nrow,ncol,figsize=(70*ncol,80*nrow))
fig.subplots_adjust(hspace=0.2, wspace=0.12) # Adjust the bottom margin, vertical space, and horizontal space
#fig.suptitle('Histograms of Hourly Statistics for each GRU', fontsize=40,y=1.0)
    
def run_loop(i,var,mx,rep,stat):
    r = i//ncol
    c = i-r*ncol
    if rep == 1: stat = 'avge'
    if rep == 2: stat = 'maxe'
    stat0 = stat
    if stat == 'rmse' or stat == 'kgem' or stat == 'avge': 
        if var == 'wallClockTime': stat0 = 'mean'
        statr = 'mean_ben'
    if stat == 'rmnz':
        if var == 'wallClockTime': stat0 = 'mnnz'
        statr = 'mnnz_ben'
    if stat == 'maxe': 
        if var == 'wallClockTime': stat0 = 'amax'
        statr = 'amax_ben'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx
    else:
        mx = 0.0
        mn = 1.0
        s_rel = summa[method_name[0]][var].sel(stat=statr)
        for m in method_name:
            s = summa[m][var].sel(stat=stat0)
            if do_rel and var != 'wallClockTime': s = s/s_rel
            if stat == 'maxe': s = np.fabs(s) # make absolute value norm
            mx = max(s.max(),mx)
            if stat == 'kgem': mn = min(s.min(),mn)

    # Data
    s_rel = summa[method_name[0]][var].sel(stat=statr)
    for m in method_name:
        s = summa[m][var].sel(stat=stat0)
        if do_rel and var != 'wallClockTime': s = s/s_rel

        if var == 'scalarTotalET' and not do_rel:
            if stat =='rmse' or stat =='rmnz' or stat=='mean': s = s*31557600 # make annual total
            if stat =='maxe': s = s*3600 # make hourly max
        if var == 'averageRoutedRunoff' and not do_rel:
            if stat =='rmse' or stat =='rmnz' or stat=='mean': s = s*31557600*1000 # make annual total
            if stat =='maxe': s = s*3600*1000 # make hourly max           
        if stat == 'maxe': s = np.fabs(s) # make absolute value norm
        range = (0,mx)
        if stat=='kgem' and var!='wallClockTime' : range = (mn,1)
        if do_hist: 
            np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=range)
        else: #cdf
            sorted_data = np.sort(np.fabs(s))
            valid_data = sorted_data[~np.isnan(sorted_data)]
            yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
            axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=2.0)
            axs[r,c].set_xlim(range)  # Replace xmin and xmax with the desired limits


    if stat0 == 'rmse': stat_word = 'RMSE'
    if stat0 == 'rmnz': stat_word = 'RMSE' # no 0s'
    if stat0 == 'maxe': stat_word = 'max abs error'
    if stat0 == 'kgem': stat_word = 'KGE"'
    if stat0 == 'mean': stat_word = 'mean'
    if stat0 == 'mnnz': stat_word = 'mean' # no 0s'
    if stat0 == 'amax': stat_word = 'max'
    if stat0 == 'avge': stat_word = 'mean abs error'
    
    if statr == 'mean_ben': statr_word = 'mean'
    if statr == 'mnnz_ben': statr_word = 'mean' # no 0s'
    if statr == 'amax_ben': statr_word = 'max'
    
    if c==0: axs[r,c].legend(plt_name)
    titl = plt_titl[i]
    if no_snow: titl = titl + ' (snow-free GRUs)'
    if rep>0: titl = titl #+ ' '+ stat_word
    axs[r,c].set_title(titl)
    if stat=='rmse' or stat=='rmnz' or stat=='maxe' or stat=='mean' or stat=='avge': axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl[i]))
    if stat=='kgem': axs[r,c].set_xlabel(stat_word)
    if do_rel and var!='wallClockTime': axs[r,c].set_xlabel('relative '+ stat_word)

    if do_hist: 
        axs[r,c].set_ylabel('GRU count')
        if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
 
    else:
        axs[r,c].set_ylabel('cumulative distribution')
        if(c>=1): axs[r, c].set_ylabel('')
        axs[r,c].set_ylim([0.0, 1.0])
        axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
        if mx<1: # Rotate x-axis labels
            axs[r, c].tick_params(axis='x', rotation=45)

def run_loopb(i,var,mx,rep,stat2):
    r = (i+len(use_vars))//ncol
    c = (i+len(use_vars))-r*ncol
    stat0 = np.copy(stat2)
    if rep == 1: stat0 = 'mean'
    if rep == 2: stat0 = 'amax'
        
    if 'zoom' in fig_fil:
        mx = mx
        mn = mx*1e-9
        if any(substring in var for substring in ['VegNrg', 'SnowNrg', 'SoilNrg']):
            mn = mx*1e-9
        if var=='wallClockTime': mn = 0.0
        if fix_units_soil and 'Soil' in var:
            mn = mn*3600*3.0 # mult by time step and depth to get storage
            mx = mx*3600*3.0
            if 'Nrg' in var:
                mn=mn*1e-3
                mx=mx*1e-3
    else:
        mx = 0.0
        mn = 1.0
        for m in method_name2:
            # Get the statistics, remove 9999 (should be nan, but just in case)
            s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
            mx = max(s.max(),mx)
            mn = min(s.min(),mn)

    # Data
    for m in method_name2:
        s = summa1[m][var].sel(stat=stat0).where(lambda x: x != 9999)
        if fix_units_soil and 'Soil' in var: 
            s = s*3600*3.0 # mult by time step and depth to get storage
            if 'Nrg' in var: s = s*1e-3

        range = (mn,mx)
        if do_hist: 
            np.fabs(s).plot.hist(ax=axs[r,c], bins=num_bins,histtype='step',zorder=0,label=m,linewidth=2.0,range=range)
        else: #cdf
            sorted_data = np.sort(np.fabs(s))
            valid_data = sorted_data[~np.isnan(sorted_data)]
            yvals = np.arange(len(valid_data)) / float(len(valid_data) - 1)
            axs[r,c].plot(valid_data, yvals, zorder=0, label=m, linewidth=2.0)
            axs[r,c].set_xlim(range)  # Replace xmin and xmax with the desired limits


    if stat0 == 'mean': stat_word = 'mean abs balance'
    if stat0 == 'amax': stat_word = 'max abs balance'

    if c==0: axs[r,c].legend(plt_name2)
    titl = plt_titl2[i]
    if no_snow: titl = titl + ' (snow-free GRUs)'
    if rep>0: titl = titl #+ ' '+ stat_word
    axs[r,c].set_title(titl)
    axs[r,c].set_xlabel(stat_word + ' [{}]'.format(leg_titl2[i]))   

    if do_hist: 
        axs[r,c].set_ylabel('GRU count')
        if(c==1): axs[r, c].set_ylabel('')
        if var != 'wallClockTime' and not run_local: axs[r,c].set_ylim([0, 25000])
 
    else:
        axs[r,c].set_ylabel('cumulative distribution')
        if(c>=1): axs[r, c].set_ylabel('')
        axs[r,c].set_ylim([0.0, 1.0])
        axs[r,c].set_xscale('log') #log x axis
        if var=='wallClockTime': 
            axs[r,c].set_xscale('function', functions=(power_transform, np.power)) #log x axis
            axs[r, c].tick_params(axis='x', rotation=45) # Rotate x-axis labels for subplot

if len(use_vars) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars,maxes,rep)): 
        run_loop(i,var,mx,rep,stat)
if len(use_vars2) > 0:
    for i,(var,mx,rep) in enumerate(zip(plot_vars2,maxes2,rep2)): 
        run_loopb(i,var,mx,rep,stat2)

# Remove the extra subplots
if (len(plot_vars)+len(plot_vars2)) < ncol*nrow:
    for i in range((len(plot_vars)+len(plot_vars2)),ncol*nrow):
        r = i//ncol
        c = i-r*ncol
        fig.delaxes(axs[r, c])

# Save
plt.savefig(viz_dir/fig_fil, bbox_inches='tight', transparent=False)
