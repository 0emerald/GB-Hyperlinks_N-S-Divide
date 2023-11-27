#!/bin/bash

# Script to create a lookup file of MSOA and region

cut -d, -f 9,10,12,13 PCD_OA_LSOA_MSOA_LAD_FEB20_UK_LU.csv > MSOA_LAD_UK_lookup.csv

# to execute, navigate in terminal to where this bash script is saved and then execute
# chmod +x MSOA_LAD_UK_lookup_bash.sh
# ./MSOA_LAD_UK_lookup_bash.sh

# not entirely sure this works, but copy and paste line 5 into terminal and run it does the job. 
