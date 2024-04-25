link_correct_files_2022_to_2023(){
    filename="$1"
    source="/data/pragallva/2022_repeat_ERA5/post_processing/parameters/$filename"
    destination="/data/pragallva/2023_repeat_ERA5/post_processing/parameters/$filename"
    ln -s "$source" "$destination"
}

# declare -a filenames=('A_N.hdf5' 'EPy1_N.hdf5' 'EPy1a_N.hdf5' 'EPy2_N.hdf5' 'EPy2a_N.hdf5' 'F1_N.hdf5' 'F2_N.hdf5' 'F3_N.hdf5' 'U_N.hdf5' 'Uref_N.hdf5' 'dEPz_dz_N.hdf5')
# for file in "${filenames[@]}"; do
#     link_correct_files_2022_to_2023 "$file"
#     echo $file
# done

link_blocking_tracked_files_2022_to_2023(){
    source="$1"       ##"/data/pragallva/2022_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/"
    destination="$2"  ##"/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/"
    ln -s "$source" "$destination"
}

# link_blocking_tracked_files_2022_to_2023 "/data/pragallva/2022_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/" "/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/"
link_blocking_tracked_files_2022_to_2023 "/data/pragallva/2022_repeat_ERA5/post_processing/full_fields/" "/data/pragallva/2023_repeat_ERA5/post_processing/"

