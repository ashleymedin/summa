#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, 6 m roots drawing on MODFLOW's EVT, without lateral flow;
# the distributed half of run_sagehen1_wet_deeproot.sh. See domain_sagehen9/README.md.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6_deeproot.config \
                      ex-gwf-sagehen-ss \
                      domain_sagehen9/settings/SUMMA/fileManager_deeproot.txt
