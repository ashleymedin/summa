#!/bin/bash
# domain_sagehen9 again, with the 9 reaches routed by mizuRoute; see domain_sagehen9/README.md.
# Needs bin/summa_modflow6_mizuroute.exe, built with -DUSE_MODFLOW6=ON -DUSE_MIZUROUTE=ON.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      -t domain_sagehen9/settings/mizuroute/mizu_control_sagehen9.toml \
                      ex-gwf-sagehen \
                      domain_sagehen9/settings/SUMMA/fileManager_mizuroute.txt \
                      ../../../bin/summa_modflow6_mizuroute.exe
