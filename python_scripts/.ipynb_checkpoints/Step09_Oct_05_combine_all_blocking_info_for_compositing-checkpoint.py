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
import scipy as sc
from scipy.interpolate import interp1d

NH =''
SH ='S'

import matplotlib
cmap1 = matplotlib.colors.LinearSegmentedColormap.from_list("", \
              [ "navy", "dodgerblue", "PowderBlue", "white", "khaki", "orange", "darkred"])


from tqdm import tqdm


####### Need to run and check this #######


def combine(lat_region = 'lat45-60', total_data_dictionary  = {}, BLOCK_DIVISIONS='DELETE_REDUNDANT_BLOCKS'):

    source = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/DJF_Ac=65/banded_data/'

    for lon_region in (['lonp120-p180', 'lonm180-m120', 'lonm120-m60', \
                        'lonm60-m0',    'lonp0-p60',    'lonp60-p120']):
        
        logging_object.write('On lon_band = %s'%(lon_region))
        ## ----> Right now I only have data for this latitude band --> This can also be put into a loop later
        data_source0         =  source+'/'+lon_region+'/'+lat_region+'/'+BLOCK_DIVISIONS+'/block_evolution_29_days/'  

        if os.path.exists(data_source0):
        
            blocking_length_days =  os.listdir(data_source0)

            total_data_dictionary[lon_region+'_'+lat_region] = {}

            for block_length in tqdm(blocking_length_days, desc=lon_region):
                data_source = data_source0+'/'+block_length+'/'
                total_data_dictionary[lon_region+'_'+lat_region][block_length] = {}

                for FIELD in os.listdir(data_source):
                    data_source2 = data_source+'/'+FIELD+'/'
                    total_data_dictionary[lon_region+'_'+lat_region][block_length][FIELD] = {}

                    if BLOCK_DIVISIONS=='DELETE_REDUNDANT_BLOCKS':
                        block_type = 'all_blocks'   
                        file_names = data_source2+'/*hkl'
                        total_data_dictionary[lon_region+'_'+lat_region][block_length][FIELD][block_type]= {}
                        
                        for file in glob.glob(file_names):
                            total_data_dictionary[lon_region+'_'+lat_region][block_length]\
                                                  [FIELD][block_type][file.split('/')[-1].split('.hkl')[0]] = hkl.load(file)
                            
                    else:
                        for block_type in os.listdir(data_source2):
                            file_names = data_source2+'/'+block_type+'/*hkl'
                            
                            total_data_dictionary[lon_region+'_'+lat_region][block_length][FIELD][block_type]= {}
                            
                            for file in glob.glob(file_names):
                                total_data_dictionary[lon_region+'_'+lat_region][block_length]\
                                                     [FIELD][block_type][file.split('/')[-1].split('.hkl')[0]] = hkl.load(file)

    return total_data_dictionary


if __name__ == "__main__":
    
    current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
    today = date.today()
    cur_date = today.strftime("%d.%m.%Y")
    logfile  = '%s-%s'%(current_file_name, cur_date)

    logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')

    for BLOCK_DIVISIONS in ['DELETE_REDUNDANT_BLOCKS', 'BlOCK_TYPES_DIVISION_LWA_based', 'BlOCK_TYPES_DIVISION_Z300_based']:
    
        start = ti.time()
        
        logging_object.write('#########################################')
        logging_object.write('--> Working on block division now --> %s'%(BLOCK_DIVISIONS))
        logging_object.write('-----------------------------------------')
        total_data_dictionary = {}
        for lat_band in ['lat20-30', 'lat30-45', 'lat45-60', 'lat60-75', 'lat75-90']:

                logging_object.write('On lat_band = %s'%(lat_band))
                logging_object.write('---    ----   ----   ---')
                total_data_dictionary = combine(lat_band, total_data_dictionary, BLOCK_DIVISIONS)

        DEST  = '/data/pragallva/2023_repeat_ERA5/post_processing/combined_data_blocking/'
        h5saveload.make_sure_path_exists(DEST)
        hkl.dump(total_data_dictionary, file_obj = DEST+'%s.hkl'%(BLOCK_DIVISIONS), track_times=True, mode='w')
        
        logging_object.write('Saved %s'%(BLOCK_DIVISIONS))
        logging_object.write('*******************************************')
        
        end = ti.time()
        time_taken = (end-start)
        logging_object.write('Time taken for this %1.3f'%(time_taken))
        logging_object.write('-----------------------------------------')
        
    logging_object.write('Great! Finished succesfully! You had a good day on %s!!'%(cur_date))
        
        