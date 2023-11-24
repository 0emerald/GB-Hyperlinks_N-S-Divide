#!/bin/bash

# Script to create a lookup file of MSOA and region

cut -d, -f 8,15 "Output_Area_to_LSOA_to_MSOA_to_Local_Authority_District_(December_2017)_Lookup_with_Area_Classifications_in_Great_Britain.csv" > MSOA_RGN11NM_lookup.csv

# to execute, navigate in terminal to where this bash script is saved and then execute
# chmod +x MSOA_REG11NM_lookup_bash.sh
# ./MSOA_REG11NM_lookup_bash.sh

