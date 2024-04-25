import numpy as np
import time as ti
import numpy.ma as ma
import hickle as hkl
import numpy as np
import sys
import logging as logg
import hickle as hkl
              
import os
import calendar
from datetime import date
import gc

import sys
sys.path.append('/data/pragallva/2023_repeat_ERA5/modules/')
import logruns as logruns
import save_and_load_hdf5_files as h5saveload



def flat_t(field):
    field_daily_linear      = np.reshape(field, (U_N.shape[0]*12*31,U_N.shape[-2],U_N.shape[-1]))
    t_indices_which_are_nan = np.argwhere(np.isnan( np.reshape(U_N, (U_N.shape[0]*12*31,U_N.shape[-2],U_N.shape[-1]) ) ))[:,0]
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
                        BLOCK_window,  lat, lon, string_dates = return_section_window( field  = FIELD,  location = blocking_details[block_event]['location'],  date = blocking_details[block_event]['date'],                                                                                                                          block_days = blocking_details[block_event]['block_days'], lon_range  = lon_range, \
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
        today    = date.today()
        cur_date = today.strftime("%d.%m.%Y")
        logfile  = '%s-%s'%(current_file_name, cur_date)
            
        logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
   

        source        = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/'
        lon_region    =  sys.argv[1] #'lonp120-p180'
        lat_region    =  sys.argv[2] #'lat45-60'

        
        CHECK_DIR  =  source+'/'+lon_region+'/'+lat_region+'/block_evolution_29_days/days4/synoptic_WA_variance_budget/'
        
        
        if os.path.exists(CHECK_DIR):   
                logging_object.write('ALREADY CALCULATED - Skipping ... \n %s'%(CHECK_DIR))
        else:    
                logging_object.write('%s \n does not exist'%CHECK_DIR)
                logging_object.write('Beginning to load all the files .....')
                route0 = '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/'
                U_N                = h5saveload.load_dict_from_hdf5(route0+'U_N.hdf5')
                lonn_N = U_N['lon']
                latt_N = U_N['lat']
                U_N    = U_N['field']

                route_filtered = '/data/pragallva/2023_repeat_ERA5/post_processing/filtered_parameters/band_pass_time_filter_004-infy_days/'
                EP_and_advective_flux_convergence_synoptic = h5saveload.load_dict_from_hdf5(route_filtered+'EPa_and_advective_flux_covergence.hdf5')
                EP_and_advective_fluxes_synoptic           = h5saveload.load_dict_from_hdf5(route_filtered+'EPa_and_advective_fluxes.hdf5')
                WA_integral_filtered                       = h5saveload.load_dict_from_hdf5(route_filtered+'wave_activity_budget_integral_a.hdf5')

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

                                        if not os.path.exists(dest+'/synoptic_integral_difference/'):

                                                ###################### Save integral differences ########################
                                                BLOCK_window_list_A      = average_over_many_blocks(FIELD = WA_integral_filtered['dA_dt'], \
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_A_diff_integral_filtered', \
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1,\
                                                                                                    blocking_details=blocking_details)[0]

                                                BLOCK_window_list_Fx     = average_over_many_blocks(FIELD = WA_integral_filtered['dF_dx'], blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_Fx_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_F1x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF1_adv_dx'], blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_F1x_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_F2x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF2_adv_dx'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_F2x_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_F3x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPx_dx'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_F3x_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_EPy    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPy_dy'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/',  name='BLOCK_window_list_EPy_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_EPz    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPz_dz'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_EPz_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                                BLOCK_window_list_residual = average_over_many_blocks(FIELD = WA_integral_filtered['residual'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral_difference/', name='BLOCK_window_list_residual_diff_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=1 )[0]

                                        else:
                                                logging_object.write('Already exists  --> \n '+dest+'/synoptic_integral_difference/')


                                        if not os.path.exists(dest+'/synoptic_integral/'):

                                                ###################### Save integral ########################
                                                BLOCK_window_list_A      = average_over_many_blocks(FIELD = WA_integral_filtered['dA_dt'], \
                                                                                                    dest=dest+'/synoptic_integral/', \
                                                                                                    name='BLOCK_window_list_A_integral_filtered', \
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0,\
                                                                                                    blocking_details=blocking_details)[0]

                                                BLOCK_window_list_Fx     = average_over_many_blocks(FIELD = WA_integral_filtered['dF_dx'], blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_Fx_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F1x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF1_adv_dx'], blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_F1x_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F2x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF2_adv_dx'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_F2x_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F3x    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPx_dx'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_F3x_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_EPy    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPy_dy'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/',  name='BLOCK_window_list_EPy_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_EPz    = average_over_many_blocks(FIELD = WA_integral_filtered['dF_EPz_dz'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_EPz_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_residual = average_over_many_blocks(FIELD = WA_integral_filtered['residual'],  blocking_details=blocking_details,\
                                                                                                    dest=dest+'/synoptic_integral/', name='BLOCK_window_list_residual_integral_filtered',\
                                                                                                    lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]


                                        else:
                                                logging_object.write('Already exists  --> \n '+dest+'/synoptic_integral/')


                                        if not os.path.exists(dest+'/synoptic_WA_budget/'):

                                                ########## Save filtered wave activity budget terms #############
                                                BLOCK_window_list_dxFfull        = average_over_many_blocks(FIELD = (EP_and_advective_flux_convergence_synoptic['dF_dx']), \
                                                                                                            dest=dest+'/synoptic_WA_budget/', \
                                                                                                            blocking_details=blocking_details,\
                                                                                                            name='BLOCK_window_list_dxF_full_synoptic',\
                                                                                                            lon_range = 60, lat_range = 15, time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF1full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF1_adv_dx'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/',name='BLOCK_window_list_dxF1_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF2full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF2_adv_dx'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/',name='BLOCK_window_list_dxF2_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF3full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPx_dx'],  blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/',name='BLOCK_window_list_dxF3_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dyEPyfull      = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPy_dy'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/',name='BLOCK_window_list_dyEPy_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dzEPzfull      = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPz_dz'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/',name='BLOCK_window_list_dzEPz_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dAdtfull       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dA_dt'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_budget/', name='BLOCK_window_list_dAdt_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]           


                                        else:
                                                logging_object.write('Already exists  --> \n '+dest+'/synoptic_WA_budget/')


                                        if not os.path.exists(dest+'/synoptic_WA_fluxes/'):

                                                ###### Save unfiltered wave activity and wave activity fluxes ######
                                                BLOCK_window_list_Afull       = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['A_N'], \
                                                                                                         dest=dest+'/synoptic_WA_fluxes/', blocking_details=blocking_details,\
                                                                                                         name='BLOCK_window_list_A_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_Ffull       = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F_adv_x']+\
                                                                                                                 EP_and_advective_fluxes_synoptic['F_EP_x'] ,  blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/',name='BLOCK_window_list_F_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F1full      = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F1_adv_x'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/',name='BLOCK_window_list_F1_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F2full      = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F2_adv_x'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/',name='BLOCK_window_list_F2_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F3full         = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F_EP_x'], blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/', name='BLOCK_window_list_F3_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F_EPy_full     = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F_EP_y'],   blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/', name='BLOCK_window_list_F_EPy_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_F_EPz_full     = average_over_many_blocks(FIELD = EP_and_advective_fluxes_synoptic['F_EP_z'],   blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_fluxes/', name='BLOCK_window_list_F_EPz_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                        else:
                                                logging_object.write('Already exists  --> \n '+dest+'/synoptic_WA_fluxes/')


                                        if not os.path.exists(dest+'/synoptic_WA_variance_budget/'):

                                                ########## Save filtered variance wave activity budget terms #############

                                                dAdt = EP_and_advective_flux_convergence_synoptic['dA_dt']

                                                BLOCK_window_list_dxFfull        = average_over_many_blocks(FIELD = (EP_and_advective_flux_convergence_synoptic['dF_dx']*dAdt), \
                                                                                                            dest=dest+'/synoptic_WA_variance_budget/', \
                                                                                                            blocking_details=blocking_details,\
                                                                                                            name='BLOCK_window_list_dxF_full_synoptic',\
                                                                                                            lon_range = 60, lat_range = 15, time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF1full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF1_adv_dx']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/',name='BLOCK_window_list_dxF1_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF2full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF2_adv_dx']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/',name='BLOCK_window_list_dxF2_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dxF3full       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPx_dx']*dAdt,  blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/',name='BLOCK_window_list_dxF3_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dyEPyfull      = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPy_dy']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/',name='BLOCK_window_list_dyEPy_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dzEPzfull      = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dF_EPz_dz']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/',name='BLOCK_window_list_dzEPz_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]

                                                BLOCK_window_list_dAdtfull       = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['dA_dt']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/', name='BLOCK_window_list_dAdtfull_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]  

                                                BLOCK_window_list_residualfull   = average_over_many_blocks(FIELD = EP_and_advective_flux_convergence_synoptic['residual']*dAdt, blocking_details=blocking_details,\
                                                                                                         dest=dest+'/synoptic_WA_variance_budget/', name='BLOCK_window_list_residual_full_synoptic',\
                                                                                                         lon_range = 60, lat_range = 15,  time_range_plus = Tp, time_range_minus = Tm, time_difference=0 )[0]  

                                                del dAdt
                                                gc.collect()

                                        ######################## ALSO SAVE VARIANCE BUDGET #########

                                        else:
                                                logging_object.write('Already exists  --> \n '+dest+'/synoptic_WA_variance_budget/')



                            else:
                                        logging_object.write('------ Lon band --> %s and Lat band --> %s --- ' %(lon_region,lat_region)) 
                                        logging_object.write('-------Day %d is not present-------' %(blocking_days))

                del EP_and_advective_flux_convergence_synoptic
                del EP_and_advective_fluxes_synoptic
                del WA_integral_filtered
                del blocking_details_day_wise
                gc.collect()

                logging_object.write('Awesome! Finished successfully')
                end0=ti.time()
                logging_object.write('<--- TOTAL Time taken --> %1.2f'%(end0-start0))
