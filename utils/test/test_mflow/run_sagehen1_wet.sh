#!/bin/bash
# One HRU coupled to MODFLOW 6 over the domain_sagehen9 event and decisions: the lumped half of
# run_sagehen9_noLatflow.sh, compared by tools/compare_lumped_distributed.py.
cd "$(dirname "$0")"
./coupler_commands.sh -c domain_sagehen1/summa_modflow6.config \
                      ex-gwf-sagehen \
                      domain_sagehen1/settings/SUMMA/fileManager_wet.txt
