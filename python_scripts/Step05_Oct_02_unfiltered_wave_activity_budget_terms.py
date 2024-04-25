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
import scipy as sc
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
import temporal_filter_for_all_years_in_one_level_with_seasonal_cycle as tf

from scipy.interpolate import interp1d
import gc

NH =''
SH ='S'

import matplotlib

def leap_day(year):
    import calendar
    if calendar.isleap(year+1979):
        return days_leap
    else:
        return days_no_leap


def time_filter_func(field, period_min=4, period_max=None, remove_seasonal_cycle=True):
    
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
        not_nan       = np.logical_not(np.isnan(field_reshape))
        not_nan_data  = field_reshape[not_nan[:,0,0],...][:,not_nan[0,:,0],:][:,:,not_nan[0,0,:]]
        tt_old        = np.arange(0, field_reshape.shape[0],1)
        tt_new        = tt_old[not_nan[:,0,0]]
        interp        = interp1d(tt_new, not_nan_data, axis=0)
        interped_data = interp(tt_old)
    else:
        interped_data = field_reshape
        
    filtered_data, freq_axis = tf.band_pass_filter_days(interped_data, seasonal_field = None, infile=infile, logging_object=None)

    if (len(field.shape) > 3):
        filtered_reshape = filtered_data.reshape(y, m, d, la, lo) ######## Reconstruct the y,m,d axis
    else:
        filtered_reshape = filtered_data ######## Reconstruction not needed
        
    return filtered_reshape



# compute_seasonal_cyle(final_field, infile = None, logging_object=None)


def compute_seasonal(field, period_min=4, period_max=None, remove_seasonal_cycle=True):
    
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
        not_nan       = np.logical_not(np.isnan(field_reshape))
        not_nan_data  = field_reshape[not_nan[:,0,0],...][:,not_nan[0,:,0],:][:,:,not_nan[0,0,:]]
        tt_old        = np.arange(0, field_reshape.shape[0],1)
        tt_new        = tt_old[not_nan[:,0,0]]
        interp        = interp1d(tt_new, not_nan_data, axis=0)
        interped_data = interp(tt_old)
    else:
        interped_data = field_reshape
        
    filtered_data, freq_axis = tf.band_pass_filter_days(compute_seasonal_cyle, infile=infile, logging_object=None)

    if (len(field.shape) > 3):
        filtered_reshape = filtered_data.reshape(y, m, d, la, lo) ######## Reconstruct the y,m,d axis
    else:
        filtered_reshape = filtered_data ######## Reconstruction not needed
        
    return filtered_reshape



def zonal_filter_func(field, N=15):
    smooth = sc.ndimage.convolve1d(field, np.ones(N)/N, axis=- 1, mode='wrap',)
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
            
        logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
        
           
        DJF = [0,1,11]  
        

        FILTER_DIREC = ""
        route     = '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/%s/'%(FILTER_DIREC)        
        dest = route
        ######################################################################################################
        ######## WAVE ACTIVITY BUDGET TERMS and THE WAVE ACTIVITY FLUXES ARE SAVED IN SMOOTHENED TERMS #######
        ######################################################################################################
        
        file_check = route+'/EP_and_advective_flux_covergence.hdf5'
        
        load_EP = True
        if not os.path.exists(file_check):
                
                
                logging_object.write('#####################################################################################')
                logging_object.write('WAVE ACTIVITY BUDGET TERMS and THE WAVE ACTIVITY FLUXES WILL BE SAVED')
                logging_object.write('-----------------------------------------------------------------------------------')

                F1_N      =  (h5saveload.load_dict_from_hdf5(route+'F1_N.hdf5')['field'])
                F2_N      =  (h5saveload.load_dict_from_hdf5(route+'F2_N.hdf5')['field'])
                F2a_N     =  (h5saveload.load_dict_from_hdf5(route+'F2a_N.hdf5')['field'])
                F3_N      =  (h5saveload.load_dict_from_hdf5(route+'F3_N.hdf5')['field'])
                F_N       =  F1_N + F2_N + F3_N
                Fa_N      =  F1_N + F2a_N + F3_N
                A_N       =  h5saveload.load_dict_from_hdf5(route+'A_N.hdf5')
                latt_N    =  A_N['lat']
                lonn_N    =  A_N['lon']
                A_N       =  A_N['field']
                
                EPy1_N    =  (h5saveload.load_dict_from_hdf5(route+'EPy1_N.hdf5')['field'])
                EPy2_N    =  (h5saveload.load_dict_from_hdf5(route+'EPy2_N.hdf5')['field'])
                EPy1a_N    =  (h5saveload.load_dict_from_hdf5(route+'EPy1a_N.hdf5')['field'])
                EPy2a_N    =  (h5saveload.load_dict_from_hdf5(route+'EPy2a_N.hdf5')['field'])
                dEPz_dz_N =  (h5saveload.load_dict_from_hdf5(route+'dEPz_dz_N.hdf5')['field'])

                dphi      = (np.diff(latt_N)).mean()

                radius         =  6378e3
                Fx_gradient    = -(np.gradient(F_N,   np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                Fax_gradient   = -(np.gradient(Fa_N,  np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                F1x_gradient   = -(np.gradient(F1_N,  np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                F2x_gradient   = -(np.gradient(F2_N,  np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                F2ax_gradient  = -(np.gradient(F2a_N, np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                F3x_gradient   = -(np.gradient(F3_N,  np.deg2rad(lonn_N), axis=-1))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None] )
                dA_dt_N_daily  =  time_tendence(A_N)

                EMF_N_divergence  =  ((EPy1_N  - EPy2_N)/(2*np.deg2rad(dphi)))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None])
                EMFa_N_divergence =  ((EPy1a_N - EPy2a_N)/(2*np.deg2rad(dphi)))/(radius*np.cos(np.deg2rad(latt_N))[None,None,None,:,None])


                heat_flux_term    =  (dEPz_dz_N)
                residual          =  (dA_dt_N_daily) - (Fx_gradient+EMF_N_divergence+heat_flux_term)
                residual_a        =  (dA_dt_N_daily) - (Fax_gradient+EMFa_N_divergence+heat_flux_term)


                EP_and_advective_flux_convergence= {'dF1_adv_dx':F1x_gradient, 'dF2_adv_dx':F2x_gradient, \
                                                    'dF_EPx_dx' :F3x_gradient, 'dF_EPy_dy' :EMF_N_divergence, \
                                                    'dF_dx': Fx_gradient, \
                                                    'dF_EPz_dz':heat_flux_term,  'dA_dt':dA_dt_N_daily, \
                                                    'residual':residual , 'lat':latt_N, 'lon':lonn_N}

                EPa_and_advective_flux_convergence= {'dF1_adv_dx':F1x_gradient, 'dF2_adv_dx':F2ax_gradient, \
                                                     'dF_EPx_dx' :F3x_gradient, 'dF_EPy_dy' :EMFa_N_divergence, \
                                                     'dF_dx': Fax_gradient, \
                                                     'dF_EPz_dz':heat_flux_term,  'dA_dt':dA_dt_N_daily, \
                                                     'residual':residual_a , 'lat':latt_N, 'lon':lonn_N}
                
                h5saveload.save_dict_to_hdf5 ( EP_and_advective_flux_convergence,  dest+'/EP_and_advective_flux_covergence.hdf5')
                h5saveload.save_dict_to_hdf5 ( EPa_and_advective_flux_convergence, dest+'/EPa_and_advective_flux_covergence.hdf5')


                logging_object.write('-----------------------------------------------------------------------------------')
                logging_object.write('EP_and_advective_flux_covergence SAVED')
                logging_object.write('#####################################################################################')
                
                load_EP = False

                #####################################################################################################
                ############################### CALCULATE FIELD INTEGRALS AND SAVE ##################################
                #####################################################################################################
                
                
        file_check = route+'/wave_activity_budget_integral.hdf5'
        load_field_integral= True
        if not os.path.exists(file_check):
            
                if load_EP:
                    EP_and_advective_flux_convergence  = h5saveload.load_dict_from_hdf5(route+'/EP_and_advective_flux_covergence.hdf5')
        
                logging_object.write('#####################################################################################')
                logging_object.write('Calculate the FIELD INTEGRALS')
                logging_object.write('-----------------------------------------------------------------------------------')

                field_integral = {}; field_integral_a = {}
                relevant_keys  = [key for key in EP_and_advective_flux_convergence.keys() \
                                       if key not in ['lat', 'lon']]
                for key in tqdm(relevant_keys):
                    field_integral[key],   orig_shape = calc_time_integral(EP_and_advective_flux_convergence[key], flatten=True)
                field_integral['original_shape']    = np.array(orig_shape)
                field_integral['lat'] = EP_and_advective_flux_convergence['lat']
                field_integral['lon'] = EP_and_advective_flux_convergence['lon']

                h5saveload.save_dict_to_hdf5 (field_integral,   route+'/wave_activity_budget_integral.hdf5')

                logging_object.write('-----------------------------------------------------------------------------------')
                logging_object.write('Saved the FIELD INTEGRALS')
                logging_object.write('#####################################################################################')
                
                ######################################################################################################
                ########################## CALCULATE CROSS CORRELATION TERMS AND SAVE ################################
                ######################################################################################################
                load_field_integral = False


                ##### Free up some memory so as not overload RAM ####
                del (field_integral)
                del (EP_and_advective_flux_convergence)
                gc.collect()


        file_check = route+'/wave_activity_budget_integral_a.hdf5'
        load_field_integral= True
        if not os.path.exists(file_check):

                if load_EP:
                    EP_and_advective_flux_convergence  = h5saveload.load_dict_from_hdf5(route+'/EPa_and_advective_flux_covergence.hdf5')
        
                logging_object.write('#####################################################################################')
                logging_object.write('Calculate the FIELD INTEGRALS')
                logging_object.write('-----------------------------------------------------------------------------------')

                field_integral = {}; field_integral_a = {}
                relevant_keys  = [key for key in EP_and_advective_flux_convergence.keys() \
                                       if key not in ['lat', 'lon']]
                for key in tqdm(relevant_keys):
                    field_integral[key],   orig_shape = calc_time_integral(EP_and_advective_flux_convergence[key], flatten=True)
                field_integral['original_shape']    = np.array(orig_shape)
                field_integral['lat'] = EP_and_advective_flux_convergence['lat']
                field_integral['lon'] = EP_and_advective_flux_convergence['lon']

                h5saveload.save_dict_to_hdf5 (field_integral,   route+'/wave_activity_budget_integral_a.hdf5')

                logging_object.write('-----------------------------------------------------------------------------------')
                logging_object.write('Saved the FIELD INTEGRALS -- correct new one.')
                logging_object.write('#####################################################################################')
                
                ######################################################################################################
                ########################## CALCULATE CROSS CORRELATION TERMS AND SAVE ################################
                ######################################################################################################
                load_field_integral = False


                ##### Free up some memory so as not overload RAM ####
                del (field_integral)
                del (EP_and_advective_flux_convergence)
                gc.collect()

        
        file_check = route+'/covar_corr/'
        if not os.path.exists(file_check):
               
                logging_object.write('###################################################################################')
                logging_object.write('Calculated the CROSS CORRELATION TERMS')
                logging_object.write('-----------------------------------------------------------------------------------')
                
                if load_field_integral:
                        field_integral   = h5saveload.load_dict_from_hdf5(route+'/wave_activity_budget_integral.hdf5')
                        field_integral_a = h5saveload.load_dict_from_hdf5(route+'/wave_activity_budget_integral_a.hdf5')

                for no_of_day_avg in range(1,5,1):
                    
                        logging_object.write('Day %d'%(no_of_day_avg))
                        correlation_dict_DJF = {}
                        covariance_dict_DJF  = {}
                        correlation_dict_DJF_a = {}
                        covariance_dict_DJF_a  = {}
                        for field_name in ['dA_dt',     'dF_dx',     'dF1_adv_dx', 'dF2_adv_dx', \
                                            'dF_EPx_dx', 'dF_EPy_dy', 'dF_EPz_dz',  'residual']:

                            logging_object.write('Calculating %s .....'%(field_name))

                            normalized_variance, covariance = cross_correlate(field_integral['dA_dt'],       \
                                                                              field_integral[field_name],    \
                                                                              no_of_day_avg = no_of_day_avg, \
                                                                              month_index   = DJF,           \
                                                                              orig_shape    = tuple(field_integral['original_shape']))

                            normalized_variance_a, covariance_a = cross_correlate(field_integral_a['dA_dt'], \
                                                                                  field_integral_a[field_name],  \
                                                                                  no_of_day_avg = no_of_day_avg, \
                                                                                  month_index   = DJF,           \
                                                                                  orig_shape    = tuple(field_integral_a['original_shape']))

                            correlation_dict_DJF.update({field_name +'_dAdt_integral': normalized_variance})
                            covariance_dict_DJF .update({field_name +'_dAdt_integral': covariance})

                            correlation_dict_DJF_a.update({field_name +'_dAdt_integral': normalized_variance_a})
                            covariance_dict_DJF_a .update({field_name +'_dAdt_integral': covariance_a})


                        correlation_dict_DJF.update({'days': no_of_day_avg, 'lat': field_integral['lat'], 'lon': field_integral['lon']})
                        covariance_dict_DJF.update( {'days': no_of_day_avg, 'lat': field_integral['lat'], 'lon': field_integral['lon']})

                        correlation_dict_DJF_a.update({'days': no_of_day_avg, 'lat': field_integral_a['lat'], 'lon': field_integral_a['lon']})
                        covariance_dict_DJF_a.update( {'days': no_of_day_avg, 'lat': field_integral_a['lat'], 'lon': field_integral_a['lon']})


                        dest1 = route+'/covar_corr/normalized_variance/'
                        h5saveload.make_sure_path_exists(dest1)
                        h5saveload.save_dict_to_hdf5(correlation_dict_DJF,   dest1+'dict_DJF_%s.hdf5'%(str.zfill('%d'%no_of_day_avg, 2)))
                        h5saveload.save_dict_to_hdf5(correlation_dict_DJF_a, dest1+'dict_DJF_%s_a.hdf5'%(str.zfill('%d'%no_of_day_avg, 2)))


                        dest2 = route+'/covar_corr/variance/'
                        h5saveload.make_sure_path_exists(dest2)
                        h5saveload.save_dict_to_hdf5(covariance_dict_DJF,    dest2+'dict_DJF_%s.hdf5'%(str.zfill('%d'%no_of_day_avg, 2)))
                        h5saveload.save_dict_to_hdf5(covariance_dict_DJF_a,  dest2+'dict_DJF_%s_a.hdf5'%(str.zfill('%d'%no_of_day_avg, 2)))        


                        logging_object.write('Saved %d in \n %s \n %s'%(no_of_day_avg, dest1, dest2))

                logging_object.write('-----------------------------------------------------------------------------------')
                logging_object.write('CROSS CORRELATION TERMS SRE CALCULATED, GOOD JOB !!!')
                logging_object.write('#####################################################################################')


                ######################################################################################################
                ########################## CALCULATE CROSS CORRELATION TERMS AND SAVE ################################
                ######################################################################################################        

        
        