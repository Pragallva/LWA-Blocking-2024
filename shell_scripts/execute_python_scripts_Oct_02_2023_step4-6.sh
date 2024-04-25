#!/bin/sh

source="/data/pragallva/2023_repeat_ERA5/codes/python_scripts/"
code1=$source/Step04_Oct_01_zsmooth_unfiltered_wave_activity_budget_terms.py
code2=$source/Step05_Oct_02_unfiltered_wave_activity_budget_terms.py
code3=$source/Step06_Oct_02_filtered_wave_activity_budget_terms.py

task(){
#    python $code1
#    python $code2
     python $code3
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