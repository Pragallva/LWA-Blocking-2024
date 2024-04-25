#!/bin/sh

##### This code aggregates daata for the block dates of 29 day length

code1=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts_Oct_2023/Step07.1_Oct_03_block_evolution_data_for_a_large_window1.py
code2=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts_Oct_2023/Step07.2_Oct_03_block_evolution_data_for_a_large_window2.py
code3=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts_Oct_2023/Step07.3_Oct_03_block_evolution_data_for_a_large_window3.py
code4=/data/pragallva/2023_repeat_ERA5/codes/python_scripts/track_blocks_scripts_Oct_2023/Step07.4_Oct_03_block_evolution_data_for_a_large_window4.py

task(){
   code=$1
   Lon_region="$2"
   Lat_region="$3"
   echo $Lon_region $Lat_region
   python $code $Lon_region $Lat_region
}


declare -a lon_region=('lonm180-m120' 'lonm120-m60' 'lonp60-p120' 'lonp120-p180' 'lonp0-p60' 'lonm60-m0')
declare -a lat_region=('lat45-60' 'lat20-30' 'lat30-45' 'lat60-75' 'lat75-90') 

for lat in "${lat_region[@]}"; do
   for lon in "${lon_region[@]}"; do
         task $code3 "$lon" "$lat" 
         task $code4 "$lon" "$lat" 
        #  task $code3 "$lon" "$lat" 
done
done

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