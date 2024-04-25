#### This script reads netcdf files from Nakamura and Huang LWA budget calculation for Northern Hemisphere and saves it in 
#### the hdf5 file format from 20N - 80N daily data from 1979-2022 and saves it in the new location.

import netCDF4 as nc
import glob
import numpy as np
import time as ti
import matplotlib.pylab as py
import numpy.ma as ma
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from datetime import datetime
import time as ti
import tarfile
import shutil

from calendar import monthrange
def number_of_days_in_month(year=2018, month=1):
    return monthrange(year, month)[1]

import sys
sys.path.append('/data/pragallva/2022_repeat_ERA5/modules/')
import logruns as logruns
import save_and_load_hdf5_files as h5saveload
import os
# os.environ["HDF5_USE_FILE_LOCKING"] = 'FALSE'
from tqdm import tqdm
import glob
import time as ti
from datetime import date

NH =''
SH ='S'

# tar_files_source='/data/nnn/ERA5/BARO_N/'
def extract_tar_make_netcdf(tar_files_source='/mnt/winds/data2/nnn/ERA5/BARO_N/', \
    year=1979, dst_path='/data/pragallva/2023_repeat_ERA5/TRASH/'):

    #### This part was written originally when the source files were tar files.
    #### For now we have netcdf data already. So no coversion needed.
    ### So dst_path == tar_files_source

    ## src_path = tar_files_source+str(year)+'.tar'
    # if os.path.isdir(dst_path) == False:
    #     os.mkdir(dst_path)

    # if not os.path.exists('%s/%d_ep2a_N.nc'%(dst_path, year)):  
    # ##### This is to check if tar is extracted or not. And extract only if is not extracted before.
    #     if src_path.endswith('.tar'):
    #         logging_object.write('*********************   %d   *********************'%(year))
    #         logging_object.write('%d tar file  extracted'%(year))
    #         tar = tarfile.open(src_path, 'r:tar')
    #         tar.extractall(dst_path)
    #         tar.close()

    # else:
    #     logging_object.write('%d netcdf exists'%(year))
        
    dst_path = tar_files_source
    return dst_path
        

def import_netcdf(file='ep2a_N', year=1979, \
                  source='/data/pragallva/2023_repeat_ERA5/TRASH/'):
    netcdf_file = source+str(year)+'_'+file+'.nc'
    
    logging_object.write('looking into %s now'%(netcdf_file))
    
    var   = nc.Dataset(netcdf_file) 
    FIELD = []
    for key in var.variables.keys():
        field = var[key][:]
        FIELD.append(field)
    FIELD = np.vstack(FIELD)
    if file != 'Urefb_N':
        FIELD = np.reshape(FIELD, (31, 4, FIELD.shape[1], FIELD.shape[2], FIELD.shape[3]))
    else:
        FIELD = np.reshape(FIELD, (31, 4, FIELD.shape[1], FIELD.shape[2], 1))
        
    latN = np.linspace(0, 90,  FIELD.shape[-2])
    lonN = np.linspace(0, 359, FIELD.shape[-1])
    
    ##### INSERT nan in those days which are supposed to not present ###
    
    m_day_lon_lat_data = np.nanmean(FIELD, axis=1).transpose((1,0,2,3))
    for month in range(1,13):
        NO_OF_DAYS = number_of_days_in_month(year, month)
        if NO_OF_DAYS < 31:
            m_day_lon_lat_data[month-1, NO_OF_DAYS:, ...] = np.nan
 
    return m_day_lon_lat_data[None,...], lonN, latN


def import_netcdf_from_tar(file='ep2a_N', year=1979, tar_files_source='/mnt/winds/data2/nnn/ERA5/BARO_N/', ):
    dst_path             = extract_tar_make_netcdf(tar_files_source, year)
    variable, lonN, latN = import_netcdf(file=file, year=year, source=dst_path)
    return variable, lonN, latN 


def combine_all_years(file='ep2a_N', start=1979, end=2022, tar_files_source='/mnt/winds/data2/nnn/ERA5/BARO_N/'):
    all_years=[]
    for year in tqdm((range(start, end+1)), desc=file):
        
        if ((year == 2022) and ('a_N' in file)):
            if file   == 'ep2a_N':
                file_new='ep2_N' 
            elif file == 'ep3a_N':
                file_new='ep3_N' 
            elif file == 'ua2a_N':
                file_new='ua2_N' 
            else:
                file_new=file
        else:
            file_new=file

        variable, lonN, latN = import_netcdf_from_tar(file_new, year, tar_files_source)
        all_years.append(variable)
    all_years = np.vstack(all_years)
    return all_years, lonN, latN


def return_lat20_80_lon_EW(source           ='/data/pragallva/2023_repeat_ERA5/', \
                           tar_files_source ='/mnt/winds/data2/nnn/ERA5/BARO_N/'):
    
    def flip_lonEW(y):
        if np.array(y).shape[-1] != 1:
            len_lon_half = np.array(y).shape[-1]//2
            z = np.append(y[...,len_lon_half:],y[...,:len_lon_half],axis=-1)
        else:
            z = y
        return z
        
    # varsi = { 'LWAb_N'  :  'A_N',      \
    #           'ua1_N'   :  'F1_N',     \
    #           'ua2_N'   :  'F2_N',     \
    #           'Ub_N'    :  'U_N',      \
    #           'Urefb_N' :  'Uref_N',   \
    #           'ep1_N'   :  'F3_N',     \
    #           'ep2_N'   :  'EPy1_N',   \
    #           'ep3_N'   :  'EPy2_N',   \
    #           'ep4_N'   :  'dEPz_dz_N',\
    #           'ep2a_N'  :  'EPy1a_N',  \
    #           'ep3a_N'  :  'EPy2a_N',  \
    #           'ua2a_N'  :  'F2a_N',    }  

    varsi = { 'aLWAb_N'  :  'aA_N',  \
              'cLWAb_N'  :  'cA_N',  }  
    
    for file in (varsi.keys()):
        
        start = ti.time()
        
        logging_object.write('%s extracting now'%(file))
        
        dest     = source+'/post_processing/parameters/'
        h5saveload.make_sure_path_exists(dest)
        filename = dest+varsi[file]+'.hdf5'
        
        if not os.path.exists(filename):
            
            all_years, lonN, latN = combine_all_years(file=file, start=1979, end=2022, \
                                    tar_files_source=tar_files_source)
            
            lonEW1     = (((lonN+180) % 360) - 180)
            lonEW      = np.append(lonEW1[len(lonEW1)//2:],lonEW1[:len(lonEW1)//2]) 
            lat_index  = np.squeeze((np.where((latN>20) & (latN<=80))))
                        
            midlat_var = flip_lonEW(all_years[...,lat_index,:])
            latN       = np.squeeze(latN[lat_index])
            dicti      = {'field':midlat_var, 'lat':latN, 'lon':lonEW} 

            dest = source+'/post_processing/parameters/' 
            h5saveload.make_sure_path_exists(dest)
            h5saveload.save_dict_to_hdf5(dicti, filename)
            logging_object.write('Saved into hdf5 -> %s'%(filename))
            
        else:
            logging_object.write('Already exists hdf5 -> %s'%(filename))
          
        end = ti.time()
        time_taken = (end-start)
        logging_object.write('Time taken for this %1.3f'%(time_taken))
        logging_object.write('############################################')
        logging_object.write('############################################')
        logging_object.write('                                            ')
        logging_object.write('                                            ')
        logging_object.write('############################################')
        logging_object.write('############################################')
            
if __name__ == "__main__":
    
    current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
    today = date.today()
    cur_date = today.strftime("%d.%m.%Y")
    logfile  = '%s-%s'%(current_file_name, cur_date)
    
    logging_object    = logruns.default_log(logfilename = logfile, log_directory = './logs/')
    return_lat20_80_lon_EW(source = '/data/pragallva/2023_repeat_ERA5/')

    # shutil.rmtree('/data/pragallva/2023_repeat_ERA5/TRASH/')
    logging_object.write('Great! Finished succesfully! You had a good day on %s!!'%(cur_date))
    
    