import pyproj
from skimage import measure
import sys
import calendar
import os
from pathlib import Path
from functools import partial
from rtree import index
import shapely.ops as ops
import shapely.geometry as shp
from scipy import stats
import matplotlib.pyplot as plt
import operator
import hickle as hkl
from scipy.stats import chisquare
from scipy.optimize import curve_fit
from scipy import optimize
import netCDF4 as nc
import glob
import numpy as np
import time as ti
import matplotlib.pylab as py
import numpy.ma as ma
import warnings
from tqdm import tqdm
from datetime import datetime
warnings.filterwarnings('ignore')


sys.path.append('/data/pragallva/2023_repeat_ERA5/modules/')
import logruns as logruns
import save_and_load_hdf5_files as h5saveload
import netcdf_utilities as ncutil
import os
os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'
from tqdm import tqdm
import glob
from PIL import Image
import copy
import itertools
from datetime import date


def return_sorted_blocks(year=1980):
        start=ti.time()
        year   = '%d.hkl'%(year)
        blocks =  hkl.load(blocking_data+year)
        end=ti.time()
        block_keys     = list(blocks.keys())
        
        start_date_object = []
        for i in range(len(block_keys)):
            y,m,d=blocks[block_keys[i]]['start_block_day']
            y =y+INITIAL_YEAR-1; m=m+1; d=d+1
            datetime_object = datetime(y,m,d)
            blocks[block_keys[i]]['start_block_day']=datetime_object
            
        new_blocks={}
        listofTuples=sorted(blocks.items(), key = lambda x: x[1]['start_block_day'])
        for elem in listofTuples :
            new_blocks.update({elem[0]: elem[1]})
               
        end=ti.time()
        print (year, 'time taken = %1.2f'%(end-start) ) #, end=', ')
        return new_blocks


def longitude_sector_and_blocking_days(coords, days, YY, dictionary={}):
    x,y = coords
    if x>180:
        x=x-360
    if (x>0    and x<=60):        
        dictionary['p0-p60'][int(days)]+=YY
        return dictionary
    
    if (x>+60   and x<=+120): 
        dictionary['p60-p120'][int(days)]+=YY
        return dictionary
    
    if (x>+120  and x<=+180):
        dictionary['p120-p180'][int(days)]+=YY
        return dictionary
    
    if (x>-180 and x<=-120):
        dictionary['m180-m120'][int(days)]+=YY        
        return dictionary
        
    if (x>-120 and x<=-60):
        dictionary['m120-m60'][int(days)]+=YY
        return dictionary
        
    if (x>-60  and x<=0):
        dictionary['m60-m0'][int(days)]+=YY
        return dictionary
    


def segregated_data_lon_bands_and_blocking_duration(year=1980, which_key='avg_no_of_blocked_days') :
    
#     vars()["blocks%i"%year] = return_sorted_blocks(year)  
    blocks=eval("blocks%i"%year)    
    
    ### Somehow the values get overwritten if copies of dictionary are defined using loops. So using this really ugly method
    band_20_30={'m180-m120':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'm120-m60': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'm60-m0':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'p0-p60':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p60-p120': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p120-p180':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]} };
    
    band_30_45={'m180-m120':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'm120-m60': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'm60-m0':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'p0-p60':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p60-p120': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p120-p180':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]} };
    
    band_45_60={'m180-m120':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'm120-m60': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'm60-m0':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'p0-p60':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p60-p120': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p120-p180':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]} };
    
    band_60_75={'m180-m120':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'm120-m60': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'm60-m0':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'p0-p60':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p60-p120': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p120-p180':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]} };
    
    band_75_90={'m180-m120':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'm120-m60': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'm60-m0':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]}, \
                'p0-p60':   {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p60-p120': {3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]},\
                'p120-p180':{3: [], 4: [], 5: [], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [], 12: [], 13: [], 14:[], 15:[], 16:[],  17:[], 18:[], 19:[], 20:[]} };    
    
    
    lat_bands_string = [   '20-30',    '30-45',    '45-60',    '60-75',    '75-90']

    
    for k in blocks.keys():
        if which_key=='start_block_date':
            necessary_fields = [blocks[k][which_key]]
        else:
            necessary_fields = blocks[k][which_key]
        for lat_bands, centroid, avg_days, necessary_field  in zip( blocks[k]['latitude_bands'], blocks[k]['envelope_centroid'], \
                                                                    blocks[k]['avg_no_of_blocked_days'], necessary_fields) :# \
            if ((avg_days<=20) and (avg_days>2)):

                  if lat_bands == '20-30':
                     band_20_30 = longitude_sector_and_blocking_days(centroid[0], avg_days, [necessary_field], band_20_30) 
                  if lat_bands == '30-45':
                     band_30_45 = longitude_sector_and_blocking_days(centroid[0], avg_days, [necessary_field], band_30_45) 
                  if lat_bands == '45-60':
                     band_45_60 = longitude_sector_and_blocking_days(centroid[0], avg_days, [necessary_field], band_45_60) 
                  if lat_bands == '60-75':
                     band_60_75 = longitude_sector_and_blocking_days(centroid[0], avg_days, [necessary_field], band_60_75) 
                  if lat_bands == '75-90':
                     band_75_90 = longitude_sector_and_blocking_days(centroid[0], avg_days, [necessary_field], band_75_90) 
    
    all_lat_bands = [band_20_30, band_30_45, band_45_60, band_60_75, band_75_90]
    
    return all_lat_bands


def return_for_all_years(func=segregated_data_lon_bands_and_blocking_duration, KEY='blocking_A'):
    
    band_20_30=[];
    band_30_45=[];
    band_45_60=[];
    band_60_75=[];
    band_75_90=[];    

    lat_bands_string = ['20-30',       '30-45',    '45-60',    '60-75',    '75-90']
    all_lat_bands    = [band_20_30, band_30_45, band_45_60, band_60_75, band_75_90]
    
    for year in range(INITIAL_YEAR,FINAL_YEAR+1):
            data=func(year, KEY)
            for l in range(len(data)):
                all_lat_bands[l].append((data[l]))
    return all_lat_bands


def RBC(block_field, lat_bands, Y, region, D, j):
    return block_field[lat_bands_string.index(lat_bands)][Y][region][D][j]

def return_blocking_details_dictionary(region   = 'm60-m0', lat_band = '45-60'):  ### Note that region denotes--> lon_region
    
    blocking_details  = {}
    blocking_details1 = {}
    total_days = range(3,20)

    M=0
    for D in total_days:
        no_of_days={}
        n=0
        for Y in range(FINAL_YEAR - INITIAL_YEAR):
            if len(all_lat_bands_loc[lat_bands_string.index(lat_band)][Y][region][D])!=0:
                for j in range(len(all_lat_bands_loc[lat_bands_string.index(lat_band)][Y][region][D])):
                    x,y = RBC(all_lat_bands_loc, lat_band, Y, region, D, j)[0]
                    if x>180:
                        x=x-360
                    block_loc=(x,y)
                    block_peakday = RBC(all_lat_bands_peakday, lat_band, Y, region, D, j)
                    block_days    = RBC(all_lat_bands_days,    lat_band, Y, region, D, j)
                    n=n+1; M=M+1
                    vars()['bloc%d'%n]={'location':block_loc, 'date':block_peakday, 'block_days':block_days}
                    vars()['bloc%d'%M]={'location':block_loc, 'date':block_peakday, 'block_days':block_days}
                    no_of_days.update({'block%d'%n:eval('bloc%d'%n)})
                    blocking_details.update({'days%d'%(D):no_of_days})
                    blocking_details1.update({'block%d'%M:eval('bloc%d'%M)})
                    
    dest = blocking_data+'/banded_data/lon%s/lat%s/'%(region,lat_band)
    os.makedirs(dest, exist_ok=True)
    hkl.dump(blocking_details,  dest+'blocking_details_day_wise.hkl', mode='w')
    hkl.dump(blocking_details1, dest+'blocking_details_all.hkl', mode='w')
                    
    return blocking_details, blocking_details1

    
if __name__ == "__main__":
        
    #### This code combines blocking data from different years and extracts/summarizes few essential info (block location, persistence and peak date)
    #### And saves the data into different lat lon bands
    
    INITIAL_YEAR = 1980
    FINAL_YEAR   = 2021
        
    blocking_data       = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/'
    lat_bands_string    = [   '20-30',    '30-45',    '45-60',    '60-75',    '75-90']
    
    for year in range(INITIAL_YEAR,FINAL_YEAR+1):
        vars()["blocks%i"%year] = return_sorted_blocks(year)      
    
    all_lat_bands_loc     = return_for_all_years(func=segregated_data_lon_bands_and_blocking_duration, KEY= 'envelope_centroid')
    all_lat_bands_peakday = return_for_all_years(func=segregated_data_lon_bands_and_blocking_duration, KEY= 'peak_block_date')
    all_lat_bands_days    = return_for_all_years(func=segregated_data_lon_bands_and_blocking_duration, KEY= 'avg_no_of_blocked_days')
    
    for lon_region in ['m180-m120', 'm120-m60','m60-m0', 'p0-p60', 'p60-p120','p120-p180']:
        for lat_region in tqdm(lat_bands_string, desc=lon_region):
            blocking_details, blocking_details1 = return_blocking_details_dictionary(region = lon_region, lat_band = lat_region)
