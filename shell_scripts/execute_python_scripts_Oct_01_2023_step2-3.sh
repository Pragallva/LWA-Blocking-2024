#!/bin/sh

code1=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/Step02_Oct_01_only_save_filtered_terms_one_by_one_seasonal.py
code2=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/Step03_Oct_01_only_save_filtered_terms_one_by_one.py

task(){
   fieldname="$1";
   python $code1 $fieldname
   python $code2 $fieldname
   echo $fieldname 
}

declare -a arr=("A_N.hdf5" "dEPz_dz_N.hdf5" "EPy1_N.hdf5" "EPy2_N.hdf5" "F1_N.hdf5" "F2_N.hdf5" "F3_N.hdf5" "U_N.hdf5" "EPy1a_N.hdf5" "EPy2a_N.hdf5" "F2a_N.hdf5") 

(
for field in "${arr[@]}"; do
    	task "$field"
done
)

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