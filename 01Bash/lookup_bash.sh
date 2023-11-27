#!/bin/bash

# Script to create a lookup file of postcode and host
# could include year aswell, but it means my laptop runs out of memory

# Uses    all.1.csv    file

cut -d, -f 2,4 all.1.csv > lookup_all.csv

# to execute, navigate in terminal to where this bash script is saved and then execute
# chmod +x lookup_bash.sh
# ./lookup_bash.sh

