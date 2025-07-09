#!/bin/bash
# Subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
#
# Inside forcingFile_out need to change forcingPath to = desForcingPath
# Inside forcingFile_out need to change initConditionFile, attributeFile, trialParamFile to = _${GRU_id} versions
# written originally by A. Van Beusekom 2025, hardwired paths run on Anvil

module load gcc/11.2.0
module load nco/4.9.3

# GRU want to subset, change to do another GRU
#just set this as GRU_nc from output error of G* on output file (here the bad file is G069381)

GRU_file=69381

GRU_id=$((GRU_file-1))
HRU_id=$GRU_id
echo "GRU order in file is ${GRU_file} starting count at 1; id is ${GRU_id} in netcdf file"

# top paths, change these to yours
homeDir=$HOME/
desforceDir=$PROJECT/users/ashley-vanbeusekom/basin_forcing/
inpforceDir=$PROJECT/users/ashley-vanbeusekom/summaNorthAmerica_forcing/
dessettingsPath=${homeDir}basin_settings/
inpsettingsPath=${homeDir}summaNorthAmerica_settings/
summa_exe=${homeDir}SummaSundials/summa/bin/summa_sundials.exe

# in paths, probably won't change, but will need to make generic fileManager_XXXXXX.txt and put in ${dessettingsPath}
# fileManager_XXXXXX.txt example provided in summa/utils/pre-processing
fileManager_in=${dessettingsPath}fileManager_XXXXXX.txt
initConditionFile_in=${inpsettingsPath}coldState.nc
attributeFile_in=${inpsettingsPath}attributes.nc
trialParamFile_in=${inpsettingsPath}trialParams.nc

# out paths, probably won't change
fileManager_out=${dessettingsPath}fileManager_${GRU_id}.txt
initConditionFile_out=${dessettingsPath}coldState_${GRU_id}.nc
attributeFile_out=${dessettingsPath}attributes_${GRU_id}.nc
trialParamFile_out=${dessettingsPath}trialParams_${GRU_id}.nc
desForcingPath=${desforceDir}basin_${GRU_id}/

# set up directory and new file Manager (will have to change things in it manually as above)
mkdir -p "$desForcingPath"
cp $fileManager_in $fileManager_out

# do the subset
ncks -d hru,$HRU_id,$HRU_id $initConditionFile_in $initConditionFile_out
echo "coldState.nc HRU ${HRU_id} subsetted"
ncks -d gru,$GRU_id,$GRU_id -d hru,$HRU_id,$HRU_id $attributeFile_in $attributeFile_out
echo "attributes.nc GRU ${GRU_id} HRU ${HRU_id} subsetted"
ncks -d hru,$HRU_id,$HRU_id $trialParamFile_in $trialParamFile_out
echo "trialParams.nc HRU ${HRU_id} subsetted"

# forcing subset has multiple files
cd $inpforceDir
for fn in NorthAmerica_remapped_*00-00-00-chunked.nc; do
    output_fn=${desForcingPath}${fn}
    ncks -d hru,$HRU_id,$HRU_id $fn $output_fn
    echo "${fn} HRU ${HRU_id} subsetted"
done
cd $homePath

# write summa command call file
runFile=${dessettingsPath}run_${GRU_id}.sh
echo "${summa_exe} -p never -s _testSumma -m ${fileManager_out} -r e" > $runFile
