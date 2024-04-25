#!/bin/sh

source="/data/pragallva/2023_repeat_ERA5/codes/python_scripts/"
code1=$source/Step08_Oct_04_remove_redundant_data_classify_into_cyclonic_anticyclonic.py
code2=$source/Step09_Oct_05_combine_all_blocking_info_for_compositing.py

task(){
#    python $code1
   python $code2
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