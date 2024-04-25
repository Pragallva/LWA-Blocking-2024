import netCDF4 as nc
import glob
import numpy as np
import time as ti
import matplotlib.pylab as py
from IPython.display import display, clear_output
import numpy.ma as ma
import warnings
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
import sys
import logging as logg

from functools import partial
import pyproj 

import logging
                
from pathlib import Path
import os
import calendar
from datetime import date
import gc

import sys
sys.path.append('/data/pragallva/2023_repeat_ERA5/modules/')
import logruns as logruns
import save_and_load_hdf5_files as h5saveload
import netcdf_utilities as ncutil



def flat_t(field):
    field_daily_linear      = np.reshape(field, (field.shape[0]*12*31,field.shape[-2],field.shape[-1]))
    t_indices_which_are_nan = np.argwhere(np.isnan( np.reshape(field, (field.shape[0]*12*31,field.shape[-2],field.shape[-1]) ) ))[:,0]
    t_indices_which_are_nan = np.sort(list(set(t_indices_which_are_nan)))
    field_daily_linear1     = np.delete(field_daily_linear, t_indices_which_are_nan, axis=0)
    return field_daily_linear1


month_names  = {0:'Jan',   1:'Feb',    2:'Mar',    3:'Apr',    4:'May',   5:'Jun',   6:'Jul',    7:'Aug',    8:'Sep',    9:'Oct',    10:'Nov',    11:'Dec'    }
days_no_leap = {  'Jan':31,  'Feb':28,   'Mar':31,   'Apr':30,   'May':31,  'Jun':30,  'Jul':31,   'Aug':31,   'Sep':30,   'Oct':31,    'Nov':30,    'Dec':31 }
days_leap    = {  'Jan':31,  'Feb':29,   'Mar':31,   'Apr':30,   'May':31,  'Jun':30,  'Jul':31,   'Aug':31,   'Sep':30,   'Oct':31,    'Nov':30,    'Dec':31 }
month_index = {'Jan' :1, 'Feb' :2,    'Mar' :3 ,'Apr' :4, 'May' :5, 'Jun' :6, 'Jul' :7, 'Aug':8,   'Sep':9, 'Oct':10,  'Nov':11,  'Dec':12}    

year_days={}

year_days_cumulative={}; month_days_leap_cumulative={}; month_days_no_leap_cumulative={}
for y in range(1979,2022):
    year_days.update({y: 366 if calendar.isleap(y) else 365})
for year in range(1979,2022):
    year_days_cumulative.update({year: np.sum([year_days[yy] for yy in range(1979,year+1)])})
year_days_cumulative.update({1978:0}) 
for mo in range(1,13):
    month_days_leap_cumulative.update({mo: np.sum([days_leap[month_names[m-1]] for m in range(1,mo+1)])})
    month_days_no_leap_cumulative.update({mo: np.sum([days_no_leap[month_names[m-1]] for m in range(1,mo+1)])})
month_days_leap_cumulative.update({0:0}); month_days_no_leap_cumulative.update({0:0})


def return_monthly_avg(field):
    field_for_that_month_that_year=np.zeros(field.shape)
    for month in range(1,13):
        for y in range(1979,2022):
            start_index             = year_days_cumulative[y-1]; end_index  = year_days_cumulative[y]
            field_for_that_year     = field[start_index:end_index, ...]
            month_days              = month_days_leap_cumulative if calendar.isleap(y) else month_days_no_leap_cumulative
            month_start_index       = month_days[month-1]; month_end_index  = month_days[month]
            field_for_that_month_that_year[y-1979,month-1,...]   =  np.nanmean(field_for_that_year[month_start_index:month_end_index, ...], axis=0)
    return field_for_that_month_that_year


def unflat_t(field):
    final_unflattened_array=np.zeros(field.shape)
    for month in range(1,13):
        for y in range(1979,2022):
            start_index                    = year_days_cumulative[y-1]; end_index  = year_days_cumulative[y]
            field_for_that_year            = field[start_index:end_index, ...]
            month_days                     = month_days_leap_cumulative if calendar.isleap(y) else month_days_no_leap_cumulative
            month_start_index              = month_days[month-1]; month_end_index  = month_days[month]
            field_for_that_month_that_year = field_for_that_year[month_start_index:month_end_index, ...]
            no_of_days_in_that_month       = month_end_index - month_start_index
            final_unflattened_array[y-1979,month-1, :no_of_days_in_that_month, ...] = field_for_that_month_that_year
            final_unflattened_array[y-1979,month-1, no_of_days_in_that_month:, ...] = np.nan            
    return final_unflattened_array


def extract_y_m_d(Y):
    def year_of(Y):
        return int(Y.split('-')[0])
    def month_of(Y):
        return int(month_index[Y.split('-')[1]])
    def day_of(Y):
        return int(Y.split('-')[2])
    
    return (year_of(Y), month_of(Y), day_of(Y))


def find_flat_index_for_the_corresponding_date(date):
    y,m,d                          = extract_y_m_d(date)
    year_start_index               = year_days_cumulative[y-1]; 
    month_days                     = month_days_leap_cumulative if calendar.isleap(y) else month_days_no_leap_cumulative
    month_start_index              = month_days[m-1]; month_end_index  = month_days[m]
    final_index=year_start_index+month_start_index+d-1
    return final_index


def find_ymd_index_for_the_corresponding_date(final_index):
        DAY = final_index
        year_key                 = min(year_days_cumulative.keys(), key=(lambda k: (DAY-year_days_cumulative[k])  if ((DAY-year_days_cumulative[k])>=0) else 1e6))
        no_of_days_into_the_year = DAY-(year_days_cumulative[year_key])
        month_days               = month_days_leap_cumulative if calendar.isleap(year_key+1) else month_days_no_leap_cumulative

        month_key                = min(month_days.keys(), key=(lambda k: (no_of_days_into_the_year-month_days[k]) if ((no_of_days_into_the_year-month_days[k])>=0) else 1e6))
        no_of_days_into_the_month = no_of_days_into_the_year-(month_days[month_key])        
        y,m,d = (year_key+1, month_key+1, no_of_days_into_the_month+1)
        return y,m,d


def ensemble_avg(X, month_index=[0,1,11]):
    return np.nanmean(np.nanmean(return_monthly_avg(X), axis=0)[month_index,...],axis=0)


def Nday_avg1(X, no_of_day_avg=1):
    if X.shape[0]<100 :
        X=flat_t(X)
    delta_X   = X[no_of_day_avg:]-  X[:-no_of_day_avg]
    delta_X   = np.append(np.zeros((no_of_day_avg, delta_X.shape[1], delta_X.shape[2])), delta_X, axis=0) 
    delta_X[0:no_of_day_avg,...]   = np.nan
    return delta_X


def string_date(y,m,d):
    return '%d-%s-%d'%(y,month_names[m-1],d)  


def return_section_window(field, location = [-16,54],  date = '2010-Feb-4', block_days=16, lon_range=10, lat_range=10, time_range=10, time_range_plus=None, time_range_minus=None):
    if field.shape[0]<100:
        field = flat_t(field)
    flat_date_index  = find_flat_index_for_the_corresponding_date(date)

    ind_lo         = np.where( (lonn_N>=(location[0]-1)) & (lonn_N<=(location[0]+1))  )
    dlon           = lonn_N[1]-lonn_N[0]
    n              = lon_range/dlon
    lon_index_list = [i%(len(lonn_N)) for i in range(int(ind_lo[0][0]-n),int(ind_lo[0][0]+n+1))]
    ind_lo         = lon_index_list

    ind_la         = np.where( (latt_N>=(location[1]-1)) & (latt_N<=(location[1]+1))  )
    dlat           = latt_N[1]-latt_N[0]
    n              = lat_range/dlat
    lat_index_list = [i for i in range(int(ind_la[0][0]-n),int(ind_la[0][0]+n+1))]
    ind_la         = lat_index_list
   
    if (time_range_plus is None) and (time_range_minus is None):
        time_range_plus=time_range; time_range_minus=time_range
    time_window      = list(np.arange(int(flat_date_index-time_range_minus)-0, int(flat_date_index+time_range_plus)+1, 1))
    
    string_dates     = [string_date( *find_ymd_index_for_the_corresponding_date(t_index) ) for t_index in time_window]
    
    A_around_block1      = np.squeeze(field[time_window,...][:,ind_la,...][...,ind_lo])
    return A_around_block1, latt_N[ind_la], lonn_N[ind_lo], string_dates



def average_over_many_blocks(FIELD, lon_range = 60, lat_range = 15, \
                             time_range_plus = 1, time_difference=1, time_range_minus=None, \
                             blocking_details='', dest='./', name='BLOCK_window_list_Fx'):
                        
            start_time = ti.time()
            
            BLOCK_window_list =[]
            date_list         =[]
            lat_list          =[]
            lon_list          =[]
            location_list     =[]

            D=time_difference;
            for block_event in blocking_details.keys():   
                print (block_event, end=',')
                
                if time_range_plus is None:
                        time_range_plus = int(blocking_details[block_event]['block_days']/2)+2
                if time_range_minus is None:
                        time_range_minus = int(blocking_details[block_event]['block_days']/2)+2

                
                if time_difference>0 :
                    
                        BLOCK_window,  lat, lon, string_dates = return_section_window(  field  = Nday_avg1(FIELD,D),  location = blocking_details[block_event]['location'],  \
                                                                                        date     = blocking_details[block_event]['date'],  block_days   = blocking_details[block_event]['block_days'], \
                                                                                        lon_range = lon_range, lat_range=lat_range,   time_range_plus  = time_range_plus,  time_range_minus = time_range_minus)
                else:
                        BLOCK_window,  lat, lon, string_dates = return_section_window(  field  = FIELD,  location = blocking_details[block_event]['location'],  date = blocking_details[block_event]['date'],                                                                                                           block_days = blocking_details[block_event]['block_days'], lon_range  = lon_range, \
                                                                                        lat_range=lat_range,time_range_plus  = time_range_plus,  time_range_minus = time_range_minus)        
                date_list.append(string_dates)
                BLOCK_window_list.append(BLOCK_window)         
                lat_list.append(lat) 
                lon_list.append(lon)
                location_list.append([blocking_details[block_event]['block_days']])
                
                del BLOCK_window
                gc.collect()
            
            dlon = lonn_N[1]-lonn_N[0]
            dlat = latt_N[1]-latt_N[0]
            lati = np.arange((0-lat_range), (0+lat_range+dlat),dlat)
            loni = np.arange((0-lon_range), (0+lon_range+dlon),dlon)
                
            #### Save here ######
            coord_info = {'lat_list':np.array(lat_list), 'lon_list':np.array(lon_list), \
                          'location_list':(location_list),'date_list':date_list, \
                          'loni':loni, 'lati':lati}
            #### Save here ######
            
            os.makedirs(dest, exist_ok=True)
            hkl.dump(np.array(BLOCK_window_list),  dest+name+'.hkl', mode='w')            
            hkl.dump(coord_info, dest+'coord_info'+'.hkl', mode='w')
                        
            #####################ß
            del BLOCK_window_list
            gc.collect()
            
            end_time = ti.time()
            time2 = (end_time - start_time)    
            logging_object.write('Succesfully calculated and saved %s data in %s' %(name,dest))
            logging_object.write('===== Time taken for this one field ======  %1.2f'%(time2))
            
            return 0,0,0,0,0,0,0
            

            
            
if __name__ == "__main__":
    
    
        current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
        today             = date.today()
        cur_date          = today.strftime("%d.%m.%Y")
        logfile           = '%s-%s'%(current_file_name, cur_date)
            
        logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
   
        logging_object.write('Beginning to load all the files .....')
        
        source        = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/'
        lon_region    =  sys.argv[1] #'lonp120-p180'
        lat_region    =  sys.argv[2] #'lat45-60'

        
        CHECK_DIR  =  source+'/'+lon_region+'/'+lat_region+'/block_evolution_29_days/days4/LWA_values_anticylconic_cyclonic/'
        
        
        if os.path.exists(CHECK_DIR):   
                logging_object.write('ALREADY CALCULATED - Skipping ... \n %s'%(CHECK_DIR))        
        
        else:
                route0 = '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/'
                U_N                = h5saveload.load_dict_from_hdf5(route0+'U_N.hdf5')
                lonn_N = U_N['lon']
                latt_N = U_N['lat']
                U_N    = U_N['field']

                LWA_values_anticylconic_cyclonic = '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/'
                LWA_cyclonic                    = h5saveload.load_dict_from_hdf5(LWA_values_anticylconic_cyclonic+'cA_N.hdf5')['field']
                LWA_anticyclonic                = h5saveload.load_dict_from_hdf5(LWA_values_anticylconic_cyclonic+'aA_N.hdf5')['field']

                logging_object.write('...Loaded all the files')


                Tp            =  14 #None   #sys.dirc[3]    ## This is the same as time_range_plus
                Tm            =  14
                blocking_details_day_wise=hkl.load(source+'/'+lon_region+'/'+lat_region+'/'+'blocking_details_day_wise.hkl')


                start0=ti.time()
                for blocking_days in range(4,20):

                            if 'days%d'%(blocking_days) in blocking_details_day_wise:

                                        logg.info('------ Lon band --> %s and Lat band --> %s --- ' %(lon_region,lat_region))
                                        logg.info('----------------Day = %d-----------------------' %(blocking_days))

                                        blocking_details=blocking_details_day_wise['days%d'%(blocking_days)]

                                        dest          = source+'/'+lon_region+'/'+lat_region+'/block_evolution_29_days/days%d/'%(blocking_days)  



                                        ###################### Save other ERA data (u, v, pv, Z) ######################

                                        if not os.path.exists(dest+'/LWA_cylonic_anticyclonic/'):


                                                    LWA_cyclonic_full       = average_over_many_blocks(FIELD = LWA_cyclonic,   blocking_details=blocking_details,\
                                                                                                       dest=dest+'/LWA_values_anticylconic_cyclonic/', name='BLOCK_window_list_LWA_cyclonic_full',\
                                                                                                       lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]


                                                    LWA_anticyclonic_full   = average_over_many_blocks(FIELD = LWA_anticyclonic,  blocking_details=blocking_details,\
                                                                                                       dest=dest+'/LWA_values_anticylconic_cyclonic/', name='BLOCK_window_list_LWA_anticyclonic_full',\
                                                                                                       lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                        else:
                                                    logging_object.write('Already exists  --> \n '+dest+'/LWA_cylonic_anticyclonic/')



                            else:
                                        logging_object.write('------ Lon band --> %s and Lat band --> %s --- ' %(lon_region,lat_region)) 
                                        logging_object.write('-------Day %d is not present-------' %(blocking_days))


                logging_object.write('Awesome! Finished successfully')
                end0=ti.time()
                logging_object.write('<--- TOTAL Time taken --> %1.2f'%(end0-start0))
