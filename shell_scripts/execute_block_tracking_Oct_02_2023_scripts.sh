#!/bin/sh

##### This code tracks the blocking events
##### And saves them based on categories (block length, lat_bands, lon_bands)

code1=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts-after_Step01/Step01.1_save_blocking_using_constant_threshold_save_extra.py
code2=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts-after_Step01/Step01.2_save_combined_blocking_info_banded_data_new.py

task(){
#    python $code1
   python $code2
}

task

# Check the exit status
# Get the full path of the script
# Extract the script name using basename
script_path=$(readlink -f "$0")
script_name=$(basename "$script_path")
if [ $? -eq 0 ]; then
    echo $script_name" ran successfully"
else
    echo $script_name" Script encountered an error"
fi