#!/bin/sh

source="/data/pragallva/2023_repeat_ERA5/codes/python_scripts/"
code=$source/Step01_Oct_01_save_netcdf_to_hdf5_files.py

task(){
   python $code
}
task

# Check the exit status
# Get the full path of the script
script_path=$(readlink -f "$0")
# Extract the script name using basename
script_name=$(basename "$script_path")
if [ $? -eq 0 ]; then
    echo $script_name" ran successfully"
else
    echo $script_name" Script encountered an error"
fi