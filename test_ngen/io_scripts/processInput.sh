python3 scripts/tile_files_by_hru.py gauge_01073000/settings/SUMMA/attributes.nc ../../../../data/gauge_01073000/forcing.nc attrib ../../../../data/gauge_01073000

python3 scripts/tile_files_by_hru.py gauge_01073000/settings/SUMMA/trialParams_default.nc ../../../../data/gauge_01073000/forcing.nc param ../../../../data/gauge_01073000

python3 scripts/tile_files_by_hru.py gauge_01073000/settings/SUMMA/coldstate.nc ../../../../data/gauge_01073000/forcing.nc init ../../../../data/gauge_01073000

#python3 scripts/tilefiless_by_hru.py gauge_01073000/forcing/SUMMA_input/summa_forcing.nc ../../../../data/gauge_01073000/forcing.nc force ../../../../data/gauge_01073000

python3 scripts/write_input.py gauge_01073000/settings/SUMMA/fileManager.txt ../../../../data/gauge_01073000/forcing.nc e gauge_01073000