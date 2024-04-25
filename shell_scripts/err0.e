nohup: ignoring input
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 451, in <module>
    save_bootstrap_raw=False, percentiles=PERCENTILES )[0]
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 271, in stats_over_many_blocks
    blocking_details, dest, name, save_bootstrap_raw=save_bootstrap_raw)
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 187, in average_over_many_blocks
    lat_range=lat_range,time_range_plus  = time_range_plus,  time_range_minus = time_range_minus)        
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 129, in return_section_window
    field = flat_t(field)
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 27, in flat_t
    field_daily_linear1     = np.delete(field_daily_linear, t_indices_which_are_nan, axis=0)
  File "<__array_function__ internals>", line 6, in delete
  File "/data/pragallva/anaconda3/lib/python3.7/site-packages/numpy/lib/function_base.py", line 4417, in delete
    new = arr[tuple(slobj)]
MemoryError: Unable to allocate 2.59 GiB for an array with shape (16071, 60, 360) and data type float64
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 441, in <module>
    save_bootstrap_raw=False, percentiles=PERCENTILES )[0]
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 271, in stats_over_many_blocks
    blocking_details, dest, name, save_bootstrap_raw=save_bootstrap_raw)
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 187, in average_over_many_blocks
    lat_range=lat_range,time_range_plus  = time_range_plus,  time_range_minus = time_range_minus)        
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 129, in return_section_window
    field = flat_t(field)
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 27, in flat_t
    field_daily_linear1     = np.delete(field_daily_linear, t_indices_which_are_nan, axis=0)
  File "<__array_function__ internals>", line 6, in delete
  File "/data/pragallva/anaconda3/lib/python3.7/site-packages/numpy/lib/function_base.py", line 4417, in delete
    new = arr[tuple(slobj)]
MemoryError: Unable to allocate 2.59 GiB for an array with shape (16071, 60, 360) and data type float64
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 361, in <module>
    EP_and_advective_fluxes_synoptic           = h5saveload.load_dict_from_hdf5(route_filtered+'EPa_and_advective_fluxes.hdf5')
  File "/data/pragallva/2023_repeat_ERA5/modules/save_and_load_hdf5_files.py", line 19, in load_dict_from_hdf5
    return recursively_load_dict_contents_from_group(h5file, '/', track)
  File "/data/pragallva/2023_repeat_ERA5/modules/save_and_load_hdf5_files.py", line 61, in recursively_load_dict_contents_from_group
    ans[key] = item[()]
  File "h5py/_objects.pyx", line 54, in h5py._objects.with_phil.wrapper
  File "h5py/_objects.pyx", line 55, in h5py._objects.with_phil.wrapper
  File "/data/pragallva/anaconda3/lib/python3.7/site-packages/h5py/_hl/dataset.py", line 562, in __getitem__
    arr = numpy.ndarray(mshape, new_dtype, order='C')
MemoryError: Unable to allocate 2.63 GiB for an array with shape (44, 12, 31, 60, 360) and data type float64
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
Traceback (most recent call last):
  File "/data/pragallva/2023_repeat_ERA5/codes/python_scripts/bootstrapping_data_Oct_2023/Step10.1_bootstrap_block_evolution_data_for_a_large_window1_percentiles.py", line 479, in <module>
    del EP_and_advective_fluxes_synoptic
NameError: name 'EP_and_advective_fluxes_synoptic' is not defined
