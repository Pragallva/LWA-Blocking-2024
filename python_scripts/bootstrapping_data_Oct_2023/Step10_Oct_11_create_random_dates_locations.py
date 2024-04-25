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
import dictionary_utility as du
import datetime as dt
import random
from tqdm import tqdm


def sub_lattice_2d(lon, lat, lon_coordinates, lat_coordinates):
    def sub_lattice_1d(field=lat, coordinate = lat_coordinates):
        sub_lattice = field[np.squeeze(np.where((field>=coordinate[0]) * (field<=coordinate[1])))]
        return sub_lattice
    sub_lon = sub_lattice_1d(lon, lon_coordinates)
    sub_lat = sub_lattice_1d(lat, lat_coordinates)
    # sub_coordinates = np.concatenate([[(x,y) for x in sub_lon] for y in sub_lat ])
    sub_coordinates = ([(x,y) for x in sub_lon for y in sub_lat ])
    return sub_coordinates

def generate_dates_for_months_and_years(years, months):
    # rn = random.randint(1, 10)
    list_of_dates = \
           [dt.datetime(year, month, day).date() 
            for year in years
            for month in months
            if 1 <= month <= 12
            for day in range(\
                random.randint(1, 10), \
                (dt.datetime(year + (month == 12), month % 12 + 1, 1) - dt.timedelta(days=1)).day + 1, \
                10+random.randint(1, 10))]

    string_dates = [(DAY.strftime("%Y-%b-%d")) for DAY in list_of_dates]
    return string_dates

generated_dates = generate_dates_for_months_and_years(np.arange(1979, 2022, dtype=int), [1, 2, 12])

def random_combination_of_loc_and_dates(    lon, lat, \
                                            lon_coordinates, lat_coordinates,\
                                            years=np.arange(1979, 2022, dtype=int), \
                                            months=[1, 2, 12], bootstrap_no=1000 ):
                                        
    generated_dates = (generate_dates_for_months_and_years(years, months))
    sub_coordinates = (sub_lattice_2d(lon, lat, lon_coordinates, lat_coordinates))

    random_dates     = random.choices(generated_dates, k=bootstrap_no)
    random_sub_coord = random.choices(sub_coordinates, k=bootstrap_no)

    coord_time_loc  = [[str_date, coord] for str_date, coord in zip(random_dates, random_sub_coord)]

    bootsrapping_coord_time_loc = {'days4':\
                                    {'block%i'%(i+1):{  'block_days':4, \
                                                        'date':coord_time_loc[i][0], \
                                                        'location': coord_time_loc[i][1]} \
                                                        for i in range(len(coord_time_loc))}}                
    return bootsrapping_coord_time_loc


if __name__ == "__main__":
    
    A_N = h5saveload.load_dict_from_hdf5('/data/pragallva/2023_repeat_ERA5/post_processing/parameters/A_N.hdf5')
    lat, lon   = A_N['lat'], A_N['lon']
    source     = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/'


    for lon_region in (os.listdir(source)) :
        for lat_region in tqdm(os.listdir(source+'/%s/'%(lon_region)), desc=lon_region):

            blocking_details_day_wise = hkl.load(source+'/'+lon_region+'/'+lat_region+'/'+'/DELETE_REDUNDANT_BLOCKS/blocking_details_day_wise.hkl')
            
            if len(blocking_details_day_wise) != 0 :
                    ### This condition makes sure that atleast one observation point exsists.
                    
                    split_regions = lon_region.split('-')[0]
                    sign = ('m', -1) if 'm' in split_regions else ('p', 1)
                    lon_coordinates = [sign[1]*int(lon_region.split('-')[i].split(sign[0])[-1]) for i in range(2)]
                    lat_coordinates = [int(lat_region.split('-')[0].split('lat')[1]), int(lat_region.split('-')[1])]

                    bootsrapping_coord_time_loc = random_combination_of_loc_and_dates(lon, lat, \
                                                    lon_coordinates, lat_coordinates,\
                                                    years=np.arange(1979, 2022, dtype=int), \
                                                    months=[1, 2, 12], bootstrap_no=500)  #### 1979-2021, DJF season
                    
                    h5saveload.make_sure_path_exists(source+'/'+lon_region+'/'+lat_region+'/'+'/BOOTSTRAPPING_BLOCKS/')
                    hkl.dump(bootsrapping_coord_time_loc, \
                            source+'/'+lon_region+'/'+lat_region+'/'+'/BOOTSTRAPPING_BLOCKS/blocking_details_day_wise.hkl', \
                            mode='w')


            





       
