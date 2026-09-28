# pre-processing folder
Helpful scripts for a variety of pre-processing purposes:
- `convert_summa_config_v2_v3.py`: convert SUMMA v2.x configuration to SUMMA v3.0.0
- `gen_coldstate.py`: create a vector cold state file for SUMMA from constant values
- `subsetGRU.sh`: subset out a NA HRU forcing, parameter, and attribute files where GRU matches HRU
- `SUMMA_merge_restarts_into_warmState.py`: combine split domain state files (with 2 dimensions, hru and gru)
- `create_lumped_to_hru_mapping.sh`: build the lumped-to-HRU mapping used to pass SUMMA runoff to mizuRoute

## Calibration observations

A calibration target reads its observations from a NetCDF file in one form: a `time` coordinate with
CF units stamping the end of the span each value covers, `time_bnds` giving that span, and one named
variable carrying its own units. These fetch observations from the web and write that form. They are
worked examples for the bundled test domains, not a general framework.

- `acquire_streamflow.py`: daily streamflow from the Water Survey of Canada or the USGS, in m3 s-1
- `acquire_grace_tws.py`: the JPL GRACE/GRACE-FO mascon anomaly for a basin or a GRU, in mm (downloading needs a NASA Earthdata login)
- `observation_netcdf.py`: the writer both share
