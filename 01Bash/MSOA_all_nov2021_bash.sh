#!/bin/bash

# Script to create a lookup file of postcode (pcds) and MSOA (msoa11)

# Uses    NSPL_NOV_2021_UK.csv    file

cut -d, -f 3,27 NSPL_NOV_2021_UK.csv > MSOA_UK_nov2021.csv

# to execute, navigate in terminal to where this bash script is saved and then execute
# chmod +x MSOA_all_nov2021_bash.sh
# ./MSOA_all_nov2021_bash.sh

