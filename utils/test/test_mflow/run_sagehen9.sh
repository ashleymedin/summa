#!/bin/bash
# One HRU per MODFLOW cell, 9 GRUs, lateral flow routed by SUMMA (groundwatr = modLatflow); see domain_sagehen9/README.md.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen9/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_sagehen9/settings/SUMMA/fileManager_latflow.txt
