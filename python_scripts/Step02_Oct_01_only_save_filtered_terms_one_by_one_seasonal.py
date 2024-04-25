#############################################################################################
#### This script calculates seasonal cycle by doing 4-24 months band pass time filtering ####
#############################################################################################
import netCDF4 as nc
import glob
import numpy as np
import time as ti
import matplotlib.pylab as py
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
# from rtree import index
from datetime import datetime
import glob
import matplotlib
import scipy as sc
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
import copy
import itertools
from datetime import date
import temporal_filter_for_all_years_in_one_level_with_seasonal_cycle as tf

from scipy.interpolate import interp1d
import json
  
NH =''
SH ='S'

import matplotlib

def leap_day(year):
    import calendar
    if calendar.isleap(year+1979):
        return days_leap
    else:
        return days_no_leap

def zonal_filter_func(field, N=15):
    smooth = sc.ndimage.convolve1d(field, np.ones(N)/N, axis= -1, mode='wrap',)
    return smooth

def flatten_t(field):
    y, m, d, la, lo = field.shape
    orig_shape      = field.shape
    
    field_reshape = field.reshape((y*m*d, la, lo))
    
    if np.isnan(field_reshape).any():
        ##### this will save some time from unnecesary interpolation
        not_nan       = np.logical_not(np.isnan(field_reshape))
        not_nan_data  = field_reshape[not_nan[:,0,0],...][:,not_nan[0,:,0],:][:,:,not_nan[0,0,:]]    
        tt_old        = np.arange(0, field_reshape.shape[0],1) 
        tt_new        = tt_old[not_nan[:,0,0]]
        interp        = interp1d(tt_new, not_nan_data, axis=0)
        interped_data = interp(tt_old)
    else:
        interped_data = field_reshape
    return interped_data, orig_shape

def orig_shape_func(field, orig_shape):
    return field.reshape(*orig_shape)

def time_tendence(field):
    flat_field, orig_shape = flatten_t(field)    
    dF_dt                  = np.gradient(flat_field, axis=0)/(24*3600)
    reshaped_data          = orig_shape_func(dF_dt, orig_shape)
    return reshaped_data  

def INTERP(field, Zsmooth=True, Nsmooth=15):
    flat_field, orig_shape = flatten_t(field)    
    reshaped_data          = orig_shape_func(flat_field, orig_shape)

    logging_object.write('DONE')

    if Zsmooth:
        reshaped_data = zonal_filter_func(reshaped_data, N=Nsmooth)
    return reshaped_data 


def time_filter_func(field, period_min=4,        period_max=None, \
                     remove_seasonal_cycle=True, Zsmooth=True):
    
    start = ti.time()
    
    infile = tf.default_infile
    infile['remove_seasonal_cycle'] = remove_seasonal_cycle
    infile['period_max']   =  period_max 
    infile['period_min']   =  period_min    
    infile['data_per_day'] =  1 
    infile['nt']           =  2**15  
    
    if (len(field.shape) > 3):
        y, m, d, la, lo = field.shape
        field_reshape = field.reshape((y*m*d, la, lo)) ######## Flatten the time axis
    else:
        field_reshape = field ######################## Time axis is already flattened
        
    if np.isnan(field_reshape).any():
        interped_data, orig_shape = flatten_t(field)
        logging_object.write('Interpolation done')
                
    else:
        logging_object.write('Did not require any interpolation')
        interped_data = field_reshape
        
    filtered_data, freq_axis = tf.band_pass_filter_days(interped_data, seasonal_field = None, infile=infile, logging_object=logging_object)
    logging_object.write('Filtering done')

    if (len(field.shape) > 3):
        filtered_reshape = filtered_data.reshape(y, m, d, la, lo) ######## Reconstruct the y,m,d axis ########
        logging_object.write('Did not require reshaping')
    else:
        filtered_reshape = filtered_data ################ Reconstruction not needed ################
        logging_object.write('Reshaping done')
             
    if Zsmooth:
        filtered_reshape = zonal_filter_func(filtered_reshape, N=15)
        logging_object.write('Smoothing by 15 grid points along longitude')
    else:
        logging_object.write('Not smoothed along longitude')
    
    end = ti.time()
    logging_object.write('Time taken --> %1.3f seconds'%(end-start))
    return filtered_reshape


# def compute_seasonal(field, period_min=4, period_max=None, remove_seasonal_cycle=True, Zsmooth=True):

######## THIS IS NOT WORKING FOR NOW. I AM USING A WORK AROUND TO CALCULATE THE SEASONAL CYCLE. I am just using a band pass filtered of 4 months to 24 months to capture the seasonal cycle
    
#     infile = tf.default_infile
#     infile['remove_seasonal_cycle'] = remove_seasonal_cycle
#     infile['period_max']   =  period_max 
#     infile['period_min']   =  period_min    
#     infile['data_per_day'] =  1 
#     infile['nt']           =  2**15  
    
#     if (len(field.shape) > 3):
#         y, m, d, la, lo = field.shape
#         field_reshape = field.reshape((y*m*d, la, lo)) ######## Flatten the time axis
#     else:
#         field_reshape = field ######################## Time axis is already flattened
        
#     if np.isnan(field_reshape).any():
#         not_nan       = np.logical_not(np.isnan(field_reshape))
#         not_nan_data  = field_reshape[not_nan[:,0,0],...][:,not_nan[0,:,0],:][:,:,not_nan[0,0,:]]
#         tt_old        = np.arange(0, field_reshape.shape[0],1)
#         tt_new        = tt_old[not_nan[:,0,0]]
#         interp        = interp1d(tt_new, not_nan_data, axis=0)
#         interped_data = interp(tt_old)
#     else:
#         interped_data = field_reshape
        
#     filtered_data = tf.compute_seasonal_cyle(interped_data, infile=infile, logging_object=logging_object)

#     if (len(field.shape) > 3):
#         filtered_reshape = filtered_data.reshape(y, m, d, la, lo) ######## Reconstruct the y,m,d axis
#     else:
#         filtered_reshape = filtered_data ######## Reconstruction not needed
        
#     if Zsmooth:
#         filtered_reshape = zonal_filter_func(filtered_reshape, N=15)
#         logging_object.write('Smoothing by 15 grid points along longitude')
#     else:
#         logging_object.write('Not smoothed along longitude')
        
#     return filtered_reshape



def calc_time_integral(field, flatten=True):
    #### If this function calculates the integral along time axis
    if flatten:
        flat_field, orig_shape = flatten_t(field) 
    integrated = sc.integrate.cumtrapz(flat_field, dx=24*3600, axis=0, initial=None)
    integrated = np.append(integrated[0,...][None,...], integrated, axis=0)
    return integrated, orig_shape



def ensemble_avg(X, month_index=[0,1,11]):
    Y = np.nanmean(\
        np.nanmean(\
        np.nanmean(X[:,month_index,...],    \
                                   axis=1), \
                                   axis=0), \
                                   axis=0)
    return Y


def cross_correlate(X,Y, X2=None, \
                    no_of_day_avg=1, \
                    month_index=[0,1,11], \
                    orig_shape=(43, 12, 31, 60, 360)):
    
    if X2 is None:
        X2=X
    delta_X   = X[no_of_day_avg:]- X[:-no_of_day_avg]
    delta_X2  = X2[no_of_day_avg:]-X2[:-no_of_day_avg]

    delta_Y   = Y[no_of_day_avg:]-Y[:-no_of_day_avg]
    delta_X   = np.append(delta_X, np.zeros((no_of_day_avg, delta_X.shape[1], delta_X.shape[2])), axis=0) 
    delta_X2  = np.append(delta_X2, np.zeros((no_of_day_avg, delta_X2.shape[1], delta_X2.shape[2])), axis=0) 

    delta_Y   = np.append(delta_Y, np.zeros((no_of_day_avg, delta_Y.shape[1], delta_Y.shape[2])), axis=0) 
    delta_X[-no_of_day_avg:,...]   = np.nan
    delta_X2[-no_of_day_avg:,...]  = np.nan
    delta_Y[-no_of_day_avg:,...]   = np.nan
    
    delta_X  = np.reshape(delta_X,  orig_shape)
    delta_Y  = np.reshape(delta_Y,  orig_shape)
    delta_X2 = np.reshape(delta_X2, orig_shape)
    
    correlation_numerator   = ensemble_avg(delta_X*delta_Y,   month_index)
    correlation_denominator = ensemble_avg(delta_X2*delta_X2, month_index)
    
    correlate = correlation_numerator/correlation_denominator
    
    return correlate, correlation_numerator

if __name__ == "__main__":
        
        current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
        today    = date.today()
        cur_date = today.strftime("%d.%m.%Y")
        logfile  = '%s-%s'%(current_file_name, cur_date)
    

        input_dicti = {'MAIN_FILE_NAME':sys.argv[1],\
                       'PERIOD_MIN':  4*31, \
                       'PERIOD_MAX': 24*31, \
                       'remove_seasonal_cycle':False, \
                       'route': '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/', \
                       'master_destination': '/data/pragallva/2023_repeat_ERA5/post_processing/filtered_parameters/', \
                       'Zsmooth':True}
        
        ### input_dicti['MAIN_FILE_NAME'] = A_N.hdf5
        
        if input_dicti['PERIOD_MIN'] is None:
            mini = str.zfill('0',3)
        else:
            mini = str.zfill(str(input_dicti['PERIOD_MIN']),3)
            
        if input_dicti['PERIOD_MAX'] is None:
            maxi = 'infy'
        else:
            maxi = str.zfill(str(input_dicti['PERIOD_MAX']),3)
#       input_dicti['time_filter_type'] = 'band_pass_time_filter_'+mini+'-'+maxi+'_days/'
        input_dicti['time_filter_type'] = 'seasonal'
        destination                     = input_dicti['master_destination']+'/'+input_dicti['time_filter_type']+'/'
        input_dicti['destination']      = destination
        
        logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
        
        if not os.path.exists(destination + input_dicti['MAIN_FILE_NAME']):
            
            route        =  input_dicti['route']
            A_N          =  h5saveload.load_dict_from_hdf5(route + input_dicti['MAIN_FILE_NAME'])
            latt_N       =  A_N['lat']
            lonn_N       =  A_N['lon']

            logging_object.write(route + input_dicti['MAIN_FILE_NAME'] + 'loaded')

            logging_object.write('Filter module called')
            
            FILTERED_A_N =  time_filter_func(A_N['field'], \
                                             period_min  = input_dicti['PERIOD_MIN'], \
                                             period_max  = input_dicti['PERIOD_MAX'], \
                                             remove_seasonal_cycle = input_dicti['remove_seasonal_cycle'], \
                                             Zsmooth = input_dicti['Zsmooth'])
            
            logging_object.write('Filtering done')

            filtered_A_N_dicti = {'field': FILTERED_A_N, 'lat':latt_N, 'lon':lonn_N}

            h5saveload.make_sure_path_exists(destination)
            h5saveload.save_dict_to_hdf5(filtered_A_N_dicti, destination + input_dicti['MAIN_FILE_NAME'])
            
            logging_object.write('SAVED --> \n'+destination +  input_dicti['MAIN_FILE_NAME'])
            
            with open(destination+'input_file_on_%s.txt'%(cur_date), 'a') as convert_file:
                convert_file.write('\n')
                convert_file.write(json.dumps(input_dicti, indent=2))
        else:
            logging_object.write('ALREADY EXISTS --> \n'+destination +  input_dicti['MAIN_FILE_NAME'])
            

