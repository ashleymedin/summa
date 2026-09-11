#!/bin/bash
# HRU cascade ordering test, coupled to MODFLOW 6 via groundwatr = modLatflow.
# Run from anywhere:   ./domain_cascade/run.sh
# Shares the ex-gwf-sagehen MODFLOW model and supplies its own coupler config with -c,
# which carries this case's HRU->cell map_file.
cd "$(dirname "$0")/.."
./coupler_commands.sh -c domain_cascade/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_cascade/settings/SUMMA/fileManager.txt
