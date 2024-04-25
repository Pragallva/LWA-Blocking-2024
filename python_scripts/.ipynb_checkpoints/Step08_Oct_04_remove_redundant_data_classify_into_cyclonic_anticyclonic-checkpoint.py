import netCDF4 as nc
import glob
import numpy as np
import time as ti
import matplotlib.pylab as py
from IPython.display import display, clear_output
import numpy.ma as ma
import warnings
warnings.filterwarnings('ignore')
from scipy import optimize
import numpy.ma as ma
from scipy.optimize import curve_fit
from scipy.stats import chisquare
import hickle as hkl
import operator
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import shapely.geometry as shp
import shapely.ops as ops
from rtree import index
from datetime import datetime
import glob
import matplotlib
import cartopy.crs as ccrs
import cartopy.util as cutil
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import shutil
import gc


from functools import partial
import pyproj 

import logging
import time as ti
                
from pathlib import Path
import os
import calendar

import sys
sys.path.append('/data/pragallva/2023_repeat_ERA5/modules/')
import logruns as logruns
import save_and_load_hdf5_files as h5saveload
import netcdf_utilities as ncutil
import os
# os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'
from tqdm import tqdm
import glob
from PIL import Image
import copy
import itertools
from datetime import date
import scipy as sc
from scipy.interpolate import interp1d

NH =''
SH ='S'

import matplotlib
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", \
              [ "navy", "dodgerblue", "PowderBlue", "white", "khaki", "orange", "darkred"])


# %%time
from tqdm import tqdm


####### Need to run and check this #######
####### This further removes any redundant block data and combine all blocking_info into a combined dictionary

def return_location_dates(blocking_details_day_wise):

    date_lists       = {}
    location_lists   = {}
    for key in blocking_details_day_wise.keys() : ##test3.keys():
        
        date_list = []
        location_list = []
        
        for block_no in blocking_details_day_wise[key].keys():
            block_date = blocking_details_day_wise[key][block_no]['date']

            date_list.append(blocking_details_day_wise[key][block_no]['date'])
            location_list.append(blocking_details_day_wise[key][block_no]['location'])

        date_lists[key]     = date_list
        location_lists[key] = location_list
        
    return date_lists, location_lists
       
    
def identify_deletion_index(date_list, location_list):

    def not_very_distant(X1,X2):    
        p0 = X1
        p1 = X2
        
        sin_lambda1 = np.sin(np.deg2rad(p0[0]))
        sin_lambda2 = np.sin(np.deg2rad(p1[0]))
        
        latitude_distance  = np.abs(p1[1]-p0[1])
        longitude_sin_distance = np.abs(sin_lambda2 - sin_lambda1)
        longitude_condition    = np.abs(2*np.cos( np.deg2rad((p0[0]+p1[0])/2) )*np.sin(np.deg2rad(20/2))) #### --> This is the formula for Sin(x1) - Sin(x2)
        
#         longitude_distance = np.abs(p1[0]-p0[0])

        RETURN=True if ((latitude_distance<4) and (longitude_sin_distance<longitude_condition)) else False
        return RETURN

    date_set     = list(set(date_list)) 
    repeat_indices = {}
    for D in (date_set):    
        idx = [i for i, j in enumerate(date_list) if j == D]
        if len(idx) > 1:
            repeat_indices[D]  = idx 

    DELETION_INDEX    = []
    DELETION_BLOCK_NO = []

    for date in repeat_indices.keys():
        for key_pairs in itertools.combinations(repeat_indices[date], 2):
            key1, key2 = (key_pairs)
            distance_boolean = not_very_distant(location_list[key1],location_list[key2])
            if distance_boolean is True:
                DELETION_INDEX   .append(key2)
                DELETION_BLOCK_NO.append('block%d'%(key2+1))
                
    return DELETION_INDEX, DELETION_BLOCK_NO
                
    
    
def copy_dir(lon_region, lat_region, DAYS,):

    dir_name  = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)
    dest_name = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)    
    
    if not os.path.exists(dest_name):    
        destination = shutil.copytree(dir_name, dest_name)    
    return dir_name, dest_name
    
    
def delete_indices_from_field_arrays(lon_region, lat_region, DAYS, delete_index):
    
    dest_name = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)  
    dir_name  = dest_name
    
    master_field_names = os.listdir(dir_name)
    for master_field_name in master_field_names:
        path = dir_name+'/%s/'%(master_field_name)
        for variable_name in os.listdir(path):
            
                file_name = path+variable_name
                variable  = hkl.load(file_name)
            
                if variable_name not in ['coord_info.hkl']:                   
                    new_variable = np.delete(variable, delete_index, axis=0)
                else:
                    new_variable = variable
                
                h5saveload.make_sure_path_exists(dest_name+'/%s/'%(master_field_name))
                hkl.dump(new_variable, dest_name+'/%s/%s'%(master_field_name, variable_name))
                
                
                
def combine(lat_region = 'lat45-60'):

    total_data_dictionary  = {}

    source = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/'

    for lon_region in tqdm(['lonp120-p180', 'lonm180-m120', 'lonm120-m60', \
                            'lonm60-m0',    'lonp0-p60',    'lonp60-p120']):
        
        # lat_region           = 'lat45-60' ## ----> Right now I only have data for this latitude band --> This can also be put into a loop later
        data_source0         =  source+'/'+lon_region+'/'+lat_region+'/DELETE_REDUNDANT_BLOCKS/block_evolution_29_days/'  

        blocking_length_days =  os.listdir(data_source0)

        total_data_dictionary[lon_region+'_'+lat_region] = {}

        for block_length in blocking_length_days:
            data_source = data_source0+'/'+block_length+'/'
            total_data_dictionary[lon_region+'_'+lat_region][block_length] = {}

            for FIELD in os.listdir(data_source):
                data_source2 = data_source+'/'+FIELD+'/'
                total_data_dictionary[lon_region+'_'+lat_region][block_length][FIELD] = {}

                for file in glob.glob(data_source2+'/*hkl'):

                    total_data_dictionary[lon_region+'_'+lat_region][block_length]\
                                         [FIELD][file.split('/')[-1].split('.hkl')[0]] = hkl.load(file)
                    
    return total_data_dictionary




###################### SOME FUNCTIONS TO CLASSIFY BLOCKS AS CYCLONIC OR ANTICYCLONIC ####################

def anom(y):
    z = y - np.nanmean(y)
    return z

def open_variables_and_classify_blocks_based_on_Z300_anomaly(lon_region, lat_region, DAYS, FIGURE=False):
    
    dir_name  = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)
    dest_name = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_Z300_based/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)     
    
    block_day_wise_file_name = hkl.load('/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl'%(lon_region, lat_region))
    
    AC_BLOCK_INFO_DICT={}; C_BLOCK_INFO_DICT={}; MIXED_BLOCK_INFO_DICT={}
    
    if DAYS not in block_day_wise_file_name.keys():
        
        print (lon_region, lat_region, DAYS, ' --> not present')
        
    else:
        
        print (lon_region, lat_region, DAYS, ' --> is  present!')
        
        block_day_wise_info  = block_day_wise_file_name[DAYS]
        
        BLOCK_INFO_ARRAY     = []
        ###### Make an array of dictionaries #########
        for block_no in ((block_day_wise_info.keys())):
            BLOCK_INFO_ARRAY.append(block_day_wise_info[block_no])

        h5saveload.make_sure_path_exists(dest_name)

        block_length         = int(DAYS.split('days')[-1])

        Z300_var_name = dir_name+'/%s/%s'%('u_v_pv_phi', 'BLOCK_window_list_Z300_full.hkl')
        Z300_field    = hkl.load(Z300_var_name)
        shape         = Z300_field.shape

        ilat          = shape[-2]//2 - 5; flat  = shape[-2]//2 + 5+1
        ilon          = shape[-1]//2 - 5; flon  = shape[-1]//2 + 5+1

        itime         = (shape[1]-1)//2 - block_length//2 - 1
        ftime         = (shape[1]-1)//2 + block_length//2 + 2


        master_field_names = os.listdir(dir_name)
        for master_field_name in master_field_names:
            path      = dir_name+'/%s/'%(master_field_name)
            dest_path = dest_name+'/%s/'%(master_field_name)

            relevant_variable_names = [name for name in os.listdir(path) if name !='coord_info.hkl' ]

            for variable_name in relevant_variable_names :

                    file_name = path+variable_name
                    variable  = hkl.load(file_name)

                    AC_VAR    = []
                    C_VAR     = [] 
                    MIXED_VAR = []

                    AC_BLOCK_INFO_ARRAY = []
                    C_BLOCK_INFO_ARRAY  = []
                    MIXED_BLOCK_INFO_ARRAY  = []

                    for Nevent in range(shape[0]):

                        Z300           = anom(Z300_field[Nevent, ...]/1e3)
                        mid_Z300       = np.nanmean(np.nanmean(Z300[:, ilat:flat, ilon:flon], axis=-1), axis=-1)
                        Z300_aggregate = np.nanmean(mid_Z300[itime:ftime])

                        if Z300_aggregate  >  1e-6:
                            typi = 'anticyclonic'
                            AC_VAR.append(variable[None,Nevent,...])
                            AC_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent])

#                             if variable_name == 'BLOCK_window_list_Z300_full.hkl':                        
#                                 print (typi, Z300_aggregate)

                        elif Z300_aggregate < -1e-6  :
                            typi = 'cyclonic'
                            C_VAR.append(variable[None, Nevent,...])
                            C_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent]) 

#                             if variable_name == 'BLOCK_window_list_Z300_full.hkl': 
#                                 print (typi, Z300_aggregate)
                        else:
                            typi = 'mixed'
                            MIXED_VAR.append(variable[None, Nevent,...])
                            MIXED_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent])

#                             if variable_name == 'BLOCK_window_list_Z300_full.hkl': 
#                                 print (typi, Z300_aggregate)

                    lrange = np.arange(82,90,1)
                    if AC_VAR:
                        typi = 'anticyclonic'
                        AC_VAR = np.vstack(AC_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/anticyclonic/')
                        hkl.dump(AC_VAR, dest_path+'/anticyclonic/%s'%(variable_name))

                        # if (variable_name == 'BLOCK_window_list_Z300_full.hkl' and FIGURE):
                        #     for length in range(AC_VAR.shape[0]):
                        #         fig, axs = py.subplots(1, ftime-itime+1,figsize=(60,4), sharex=True, sharey=True)
                        #         ii=-1
                        #         for T in range(itime, ftime+1):
                        #             ii=ii+1
                        #             im = axs[ii].contourf(loni, lati, AC_VAR[length, T, ...]/1e3, lrange, cmap=py.cm.Spectral_r, extend='both'); 
                        #             axs[ii].tick_params(labelsize=30)
                        #         fig.suptitle('ANTICYCLONIC', fontsize=30)
                        #         fig.colorbar(im, ax=axs.ravel())

                    if C_VAR:
                        typi   = 'cyclonic'
                        C_VAR  = np.vstack(C_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/cyclonic/')
                        hkl.dump(C_VAR, dest_path+'/cyclonic/%s'%(variable_name))

                        # if (variable_name == 'BLOCK_window_list_Z300_full.hkl' and FIGURE):
                        #     for length in range(C_VAR.shape[0]):
                        #         fig, axs = py.subplots(1, ftime-itime+1,figsize=(60,4), sharex=True, sharey=True)
                        #         ii=-1
                        #         for T in range(itime, ftime+1):
                        #             ii=ii+1
                        #             im = axs[ii].contourf(loni, lati, C_VAR[length, T, ...]/1e3, lrange, cmap=py.cm.Spectral_r, extend='both'); 
                        #             axs[ii].tick_params(labelsize=30)
                        #         fig.suptitle('CYCLONIC', fontsize=30)
                        #         fig.colorbar(im, ax=axs.ravel())

                    if MIXED_VAR:
                        typi       = 'mixed'
                        MIXED_VAR  =  np.vstack(MIXED_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/mixed/')
                        hkl.dump(MIXED_VAR, dest_path+'/mixed/%s'%(variable_name))


        if (len (MIXED_VAR)> 0):
            print ('MIXED var shape --> ', MIXED_VAR.shape) 
        if (len (C_VAR)> 0):
            print ('C var shape --> ',     C_VAR.shape)  
        if (len (AC_VAR)>0) :
            print ('AC var shape --> ',    AC_VAR.shape)  

        master_field_name = 'BLOCK_INFO'
        dest_path =  dest_name+'/%s/'%(master_field_name)

        if AC_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/anticyclonic/')        
            for block_no in range(len(AC_BLOCK_INFO_ARRAY)):
                AC_BLOCK_INFO_DICT['block%d'%(block_no+1)]   = AC_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(AC_BLOCK_INFO_DICT, dest_path+'/anticyclonic/%s'%('date_location_block_no'))


        if C_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/cyclonic/')   
            for block_no in range(len(C_BLOCK_INFO_ARRAY)):
                C_BLOCK_INFO_DICT['block%d'%(block_no+1)]     =   C_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(C_BLOCK_INFO_DICT, dest_path+'/cyclonic/%s'%('date_location_block_no'))

        if MIXED_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/mixed/')   
            for block_no in range(len(MIXED_BLOCK_INFO_ARRAY)):
                MIXED_BLOCK_INFO_DICT['block%d'%(block_no+1)] =   MIXED_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(MIXED_BLOCK_INFO_DICT, dest_path+'/mixed/%s'%('date_location_block_no'))
    
    
    return AC_BLOCK_INFO_DICT, C_BLOCK_INFO_DICT, MIXED_BLOCK_INFO_DICT
            



def open_variables_and_classify_blocks_based_on_LWA_anomaly(lon_region, lat_region, DAYS, FIGURE=False):
    
    dir_name  = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)
    dest_name = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_LWA_based/block_evolution_29_days/%s/'%(lon_region, lat_region, DAYS)     
    
    block_day_wise_file_name = hkl.load('/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl'%(lon_region, lat_region))
    
    AC_BLOCK_INFO_DICT={}; C_BLOCK_INFO_DICT={}; MIXED_BLOCK_INFO_DICT={}
    
    if DAYS not in block_day_wise_file_name.keys():
        
        print (lon_region, lat_region, DAYS, ' --> not present')
        
    else:
        
        print (lon_region, lat_region, DAYS, ' --> is  present!')
        
        block_day_wise_info  = block_day_wise_file_name[DAYS]
        
        BLOCK_INFO_ARRAY     = []
        ###### Make an array of dictionaries #########
        for block_no in ((block_day_wise_info.keys())):
            BLOCK_INFO_ARRAY.append(block_day_wise_info[block_no])

        h5saveload.make_sure_path_exists(dest_name)

        block_length         = int(DAYS.split('days')[-1])

        aLWA_var_name = dir_name+'/%s/%s'%('LWA_values_anticylconic_cyclonic', 'BLOCK_window_list_LWA_anticyclonic_full.hkl')
        cLWA_var_name = dir_name+'/%s/%s'%('LWA_values_anticylconic_cyclonic', 'BLOCK_window_list_LWA_cyclonic_full.hkl')
        aLWA_field    = hkl.load(aLWA_var_name)
        cLWA_field    = hkl.load(cLWA_var_name)
        shape         = aLWA_field.shape
        
        ilat          = shape[-2]//2 - 5; flat  = shape[-2]//2 + 5+1
        ilon          = shape[-1]//2 - 5; flon  = shape[-1]//2 + 5+1
        itime         = (shape[1]-1)//2 - block_length//2 - 1
        ftime         = (shape[1]-1)//2 + block_length//2 + 2


        master_field_names = os.listdir(dir_name)
        for master_field_name in master_field_names:
            path      = dir_name+'/%s/'%(master_field_name)
            dest_path = dest_name+'/%s/'%(master_field_name)

            relevant_variable_names = [name for name in os.listdir(path) if name !='coord_info.hkl' ]

            for variable_name in relevant_variable_names :

                    file_name = path+variable_name
                    variable  = hkl.load(file_name)

                    AC_VAR    = []
                    C_VAR     = [] 
                    MIXED_VAR = []

                    AC_BLOCK_INFO_ARRAY = []
                    C_BLOCK_INFO_ARRAY  = []
                    MIXED_BLOCK_INFO_ARRAY  = []

                    for Nevent in range(shape[0]):


                        aLWA           = anom(aLWA_field[Nevent, ...])
                        cLWA           = anom(cLWA_field[Nevent, ...])
                        mid_aLWA       = np.nanmean(np.nanmean(aLWA[:, ilat:flat, ilon:flon], axis=-1), axis=-1)
                        mid_cLWA       = np.nanmean(np.nanmean(cLWA[:, ilat:flat, ilon:flon], axis=-1), axis=-1)
                        aLWA_aggregate = np.nanmean(mid_aLWA[itime:ftime])
                        cLWA_aggregate = np.nanmean(mid_cLWA[itime:ftime])

                        LWA_percent       = 100*(aLWA_aggregate-cLWA_aggregate)/((aLWA_aggregate+cLWA_aggregate)/2.0)
                        threshold_percent = 0.01

                        # Z300           = anom(Z300_field[Nevent, ...]/1e3)
                        # mid_Z300       = np.nanmean(np.nanmean(Z300[:, ilat:flat, ilon:flon], axis=-1), axis=-1)
                        # Z300_aggregate = np.nanmean(mid_Z300[itime:ftime])

                        if LWA_percent  >  threshold_percent:
                            typi = 'anticyclonic'
                            AC_VAR.append(variable[None,Nevent,...])
                            AC_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent])


                        elif LWA_percent < -threshold_percent  :
                            typi = 'cyclonic'
                            C_VAR.append(variable[None, Nevent,...])
                            C_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent]) 

                        else:
                            typi = 'mixed'
                            MIXED_VAR.append(variable[None, Nevent,...])
                            MIXED_BLOCK_INFO_ARRAY.append(BLOCK_INFO_ARRAY[Nevent])


                    lrange = np.arange(82,90,1)
                    if AC_VAR:
                        typi = 'anticyclonic'
                        AC_VAR = np.vstack(AC_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/anticyclonic/')
                        hkl.dump(AC_VAR, dest_path+'/anticyclonic/%s'%(variable_name))

                        # if (variable_name == 'BLOCK_window_list_LWA_anticyclonic_full.hkl' and FIGURE):
                        #     for length in range(AC_VAR.shape[0]):
                        #         fig, axs = py.subplots(1, ftime-itime+1,figsize=(60,4), sharex=True, sharey=True)
                        #         ii=-1
                        #         for T in range(itime, ftime+1):
                        #             ii=ii+1
                        #             im = axs[ii].contourf(loni, lati, AC_VAR[length, T, ...], lrange, cmap=py.cm.Spectral_r, extend='both'); 
                        #             axs[ii].tick_params(labelsize=30)
                        #         fig.suptitle('ANTICYCLONIC LWA', fontsize=30)
                        #         fig.colorbar(im, ax=axs.ravel())

                    if C_VAR:
                        typi   = 'cyclonic'
                        C_VAR  = np.vstack(C_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/cyclonic/')
                        hkl.dump(C_VAR, dest_path+'/cyclonic/%s'%(variable_name))

                        # if (variable_name == 'BLOCK_window_list_LWA_cyclonic_full.hkl' and FIGURE):
                        #     for length in range(C_VAR.shape[0]):
                        #         fig, axs = py.subplots(1, ftime-itime+1,figsize=(60,4), sharex=True, sharey=True)
                        #         ii=-1
                        #         for T in range(itime, ftime+1):
                        #             ii=ii+1
                        #             im = axs[ii].contourf(loni, lati, C_VAR[length, T, ...]/1e3, lrange, cmap=py.cm.Spectral_r, extend='both'); 
                        #             axs[ii].tick_params(labelsize=30)
                        #         fig.suptitle('CYCLONIC', fontsize=30)
                        #         fig.colorbar(im, ax=axs.ravel())

                    if MIXED_VAR:
                        typi       = 'mixed'
                        MIXED_VAR  =  np.vstack(MIXED_VAR)

                        h5saveload.make_sure_path_exists(dest_path+'/mixed/')
                        hkl.dump(MIXED_VAR, dest_path+'/mixed/%s'%(variable_name))


        if (len (MIXED_VAR)> 0):
            print ('MIXED var shape --> ', MIXED_VAR.shape) 
        if (len (C_VAR)> 0):
            print ('C var shape --> ', C_VAR.shape)  
        if (len (AC_VAR)>0) :
            print ('AC var shape --> ', AC_VAR.shape)  

        master_field_name = 'BLOCK_INFO'
        dest_path =  dest_name+'/%s/'%(master_field_name)

        if AC_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/anticyclonic/')        
            for block_no in range(len(AC_BLOCK_INFO_ARRAY)):
                AC_BLOCK_INFO_DICT['block%d'%(block_no+1)]   = AC_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(AC_BLOCK_INFO_DICT, dest_path+'/anticyclonic/%s'%('date_location_block_no'))


        if C_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/cyclonic/')   
            for block_no in range(len(C_BLOCK_INFO_ARRAY)):
                C_BLOCK_INFO_DICT['block%d'%(block_no+1)]     =   C_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(C_BLOCK_INFO_DICT, dest_path+'/cyclonic/%s'%('date_location_block_no'))

        if MIXED_BLOCK_INFO_ARRAY:
            h5saveload.make_sure_path_exists(dest_path+'/mixed/')   
            for block_no in range(len(MIXED_BLOCK_INFO_ARRAY)):
                MIXED_BLOCK_INFO_DICT['block%d'%(block_no+1)] =   MIXED_BLOCK_INFO_ARRAY[block_no]
            hkl.dump(MIXED_BLOCK_INFO_DICT, dest_path+'/mixed/%s'%('date_location_block_no'))
    
    
    return AC_BLOCK_INFO_DICT, C_BLOCK_INFO_DICT, MIXED_BLOCK_INFO_DICT





if __name__ == "__main__":

    current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
    today    = date.today()
    cur_date = today.strftime("%d.%m.%Y")
    logfile  = '%s-%s'%(current_file_name, cur_date)
        
    logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
    

    ########## REMOVE REDUNDANT BLOCKS #######
    source = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/'

    logging_object.write('#####################################################################################')
    logging_object.write('DELETING REDUNDANT DATA')
    logging_object.write('-----------------------------------------------------------------------------------')

    for lon_region in (os.listdir(source)):
#         for lat_region in ['lat45-60']:
        for lat_region in ['lat45-60', 'lat20-30', 'lat30-45', 'lat60-75', 'lat75-90']:

            logging_object.write('Delete redundant data --> %s - %s'%(lon_region, lat_region))

            blocking_details_day_wise = hkl.load('/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/blocking_details_day_wise.hkl'%(lon_region, lat_region))
            if 'days3' in blocking_details_day_wise.keys():
                blocking_details_day_wise.pop('days3') #### Because rest all scripts have saved only for days>=4. Those are only identified as blocks

            date_lists, location_lists = return_location_dates(blocking_details_day_wise)
            DELETION_INDICES   = {}
            DELETION_BLOCK_NOS = {}

            for key in date_lists.keys():
                DELETION_INDICES[key], DELETION_BLOCK_NOS[key] = identify_deletion_index(date_lists[key], location_lists[key])


            for DAYS in DELETION_INDICES.keys():
                copy_dir(lon_region, lat_region, DAYS,)

                if len(DELETION_INDICES[DAYS]) > 0 :
                    for deletion_elements in DELETION_BLOCK_NOS[DAYS]:
                        blocking_details_day_wise[DAYS].pop(deletion_elements)
                        print ('deleted '+lon_region+' - '+lat_region+'--> '+DAYS+', '+deletion_elements)    
                    delete_indices_from_field_arrays(lon_region, lat_region, DAYS, DELETION_INDICES[DAYS])                
                    
            h5saveload.make_sure_path_exists(\
                      '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/'%(lon_region, lat_region))
            hkl.dump(blocking_details_day_wise, \
                     '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl'%(lon_region, lat_region))

            logging_object.write('New ones saved in %s'%('DELETE_REDUNDANT_BLOCKS'))
               
    logging_object.write('Redundant data deleted')
    logging_object.write('-----------------------------------------------------------------------------------')
    


    logging_object.write('#####################################################################################')
    logging_object.write('CLASSIFY BLOCKS BASED ON Z300 ANOMALY')
    logging_object.write('-----------------------------------------------------------------------------------')

    ########### CLASSIFY BLOCKS INTO ANTICYCLONES and CYCLONES based on geopotential height anomaly #######
    source = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/'
    
    for lon_region in (os.listdir(source)):
#         for lat_region in ['lat45-60']:
        for lat_region in ['lat45-60', 'lat20-30', 'lat30-45', 'lat60-75', 'lat75-90']:
            
            logging_object.write('Classify blocks on z300 anomaly --> %s - %s'%(lon_region, lat_region))

            blocking_details_day_wise = hkl.load(\
            '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl'%(lon_region, lat_region))

            if 'days3' in blocking_details_day_wise.keys():
                blocking_details_day_wise.pop('days3') #### Because rest all scripts have saved only for days>=4. Those are only identified as blocks

            blocking_details_day_wise_division_wise = {}
            
            AC_BLOCK_INFO_DICT_WITH_DAYS={}; 
            C_BLOCK_INFO_DICT_WITH_DAYS ={}; 
            MIXED_BLOCK_INFO_DICT_WITH_DAYS={}

            for days in blocking_details_day_wise.keys():
                AC_BLOCK_INFO_DICT, C_BLOCK_INFO_DICT, MIXED_BLOCK_INFO_DICT = open_variables_and_classify_blocks_based_on_Z300_anomaly(lon_region, lat_region, days)
                
                if bool(AC_BLOCK_INFO_DICT):
                    AC_BLOCK_INFO_DICT_WITH_DAYS[days]    =  AC_BLOCK_INFO_DICT
                    
                if bool(C_BLOCK_INFO_DICT):
                    C_BLOCK_INFO_DICT_WITH_DAYS[days]     =   C_BLOCK_INFO_DICT
                                        
                if bool(MIXED_BLOCK_INFO_DICT):
                    MIXED_BLOCK_INFO_DICT_WITH_DAYS[days] = MIXED_BLOCK_INFO_DICT
                    
            blocking_details_day_wise_division_wise['anticyclonic'] =  AC_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['cyclonic']     =  C_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['mixed']        =  MIXED_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['all']          =  blocking_details_day_wise
            

            h5saveload.make_sure_path_exists(\
                      '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_Z300_based/'%(lon_region, lat_region))

            hkl.dump(blocking_details_day_wise_division_wise, \
            '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_Z300_based/blocking_details_day_wise.hkl'%(lon_region, lat_region))
    
            logging_object.write('New ones saved in %s'%('BlOCK_TYPES_DIVISION_Z300_based'))

    del (blocking_details_day_wise_division_wise)
    del (AC_BLOCK_INFO_DICT_WITH_DAYS)
    del (C_BLOCK_INFO_DICT_WITH_DAYS)
    del (MIXED_BLOCK_INFO_DICT_WITH_DAYS)
    del (AC_BLOCK_INFO_DICT)
    del (C_BLOCK_INFO_DICT)
    del (MIXED_BLOCK_INFO_DICT)
    gc.collect()
    
    logging_object.write('Blocks classified based on Z300 anomaly')
    logging_object.write('-----------------------------------------------------------------------------------')



    logging_object.write('#####################################################################################')
    logging_object.write('CLASSIFY BLOCKS BASEN ON LWA ANOMALY')
    logging_object.write('-----------------------------------------------------------------------------------')

    ########### CLASSIFY BLOCKS INTO ANTICYCLONES and CYCLONES based on LWA anomaly #######
    source = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/'
    
    for lon_region in (os.listdir(source)):
#         for lat_region in ['lat45-60']:
        for lat_region in ['lat45-60', 'lat20-30', 'lat30-45', 'lat60-75', 'lat75-90']:
            
            logging_object.write('Classify blocks on LWA anomaly --> %s - %s'%(lon_region, lat_region))

            blocking_details_day_wise = hkl.load(\
            '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl'%(lon_region, lat_region))

            if 'days3' in blocking_details_day_wise.keys():
                blocking_details_day_wise.pop('days3') #### Because rest all scripts have saved only for days>=4. Those are only identified as blocks

            blocking_details_day_wise_division_wise = {}
            
            AC_BLOCK_INFO_DICT_WITH_DAYS={}; 
            C_BLOCK_INFO_DICT_WITH_DAYS ={}; 
            MIXED_BLOCK_INFO_DICT_WITH_DAYS={}

            
            for days in blocking_details_day_wise.keys():
                AC_BLOCK_INFO_DICT, C_BLOCK_INFO_DICT, MIXED_BLOCK_INFO_DICT = open_variables_and_classify_blocks_based_on_LWA_anomaly(lon_region, lat_region, days)
                
                if bool(AC_BLOCK_INFO_DICT):
                    AC_BLOCK_INFO_DICT_WITH_DAYS[days]    =   AC_BLOCK_INFO_DICT
                    
                if bool(C_BLOCK_INFO_DICT):
                    C_BLOCK_INFO_DICT_WITH_DAYS[days]     =   C_BLOCK_INFO_DICT
                                        
                if bool(MIXED_BLOCK_INFO_DICT):
                    MIXED_BLOCK_INFO_DICT_WITH_DAYS[days] =   MIXED_BLOCK_INFO_DICT
                    
            blocking_details_day_wise_division_wise['anticyclonic'] =  AC_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['cyclonic']     =  C_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['mixed']        =  MIXED_BLOCK_INFO_DICT_WITH_DAYS
            blocking_details_day_wise_division_wise['all']          =  blocking_details_day_wise
            

            h5saveload.make_sure_path_exists(\
                      '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_LWA_based/'%(lon_region, lat_region))

            hkl.dump(blocking_details_day_wise_division_wise, \
            '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable//DJF_Ac=65/banded_data/%s/%s/BlOCK_TYPES_DIVISION_LWA_based/blocking_details_day_wise.hkl'%(lon_region, lat_region))
 
            logging_object.write('New ones saved in %s'%('BlOCK_TYPES_DIVISION_LWA_based'))

    del (blocking_details_day_wise_division_wise)
    del (AC_BLOCK_INFO_DICT_WITH_DAYS)
    del (C_BLOCK_INFO_DICT_WITH_DAYS)
    del (MIXED_BLOCK_INFO_DICT_WITH_DAYS)
    del (AC_BLOCK_INFO_DICT)
    del (C_BLOCK_INFO_DICT)
    del (MIXED_BLOCK_INFO_DICT)
    gc.collect()

    print ('Blocks classified .......')

    logging_object.write('Blocks classified based on LWA anomaly')
    logging_object.write('-----------------------------------------------------------------------------------')

logging_object.write('Great! Finished succesfully! You had a good day on %s!!'%(cur_date))

#     for lat_band in ['lat45-60']:
#         SOURCE   = '/data/pragallva/2022_repeat_ERA5/post_processing/combined_data_blocking_BLOCKING_DIVISION/'
#         h5saveload.make_sure_path_exists(SOURCE)
       
#         total_data_dictionary = combine(lat_region)
        
#         hkl.dump(total_data_dictionary, file_obj = SOURCE+'%s.hdf5'%(lat_band), track_times=True, mode='w')