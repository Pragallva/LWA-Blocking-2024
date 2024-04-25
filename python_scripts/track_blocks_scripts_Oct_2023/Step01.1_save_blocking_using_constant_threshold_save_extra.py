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
# import shapely.speedups as speedup
# speedup.disable()
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


def latavg(x, year=0, month=0, day=0):
    x0 = x[year, month, day, ...]
    x1 = x0[::2, :]
    x2 = x0[1::2, :]
    x = (x1+x2)*0.5
    return x[:-1]


def science_blocks():
    py.plot(-9, 45,   'r*',   markersize=15,)  # mfc='none')
    py.plot(-9, 45,   'wo',   markersize=2, )  # mfc='none')
    py.plot(-147, 42, 'r*',   markersize=15,)  # mfc='none')
    py.plot(-147, 42, 'wo',   markersize=2, )  # mfc='none')


def area(x, y):
    a = 0.5*np.sum(y*np.gradient(x) - x*np.gradient(y))
    a = np.abs(a)
    return a


def centroid(x, y):
    x1 = np.sum(x)/len(x)
    y1 = np.sum(y)/len(y)
    return x1, y1


def length(x, y):
    dx = np.max(x)-np.min(x)
    dy = np.max(y)-np.min(y)
    return dx, dy


def add_cyclic(YY, x=None):
    x=lon
    YY = np.append(YY, YY[..., :int(len(x)/2)], axis=-1)
    return YY


def avg_lat_lon_extent(x, y, centers):
    xcoords = np.copy(x)
    ycoords = np.copy(y)
    x_max = np.max(np.abs(xcoords-np.mean(xcoords)))
    y_max = np.max(np.abs(ycoords-np.mean(ycoords)))
    return x_max, y_max


# logging_object.write('Defined many functions')
# logging_object.write('Defined leap year, month, days')


def leap_day(year):
    if calendar.isleap(year+ZEROTH_YEAR):
        return days_leap
    else:
        return days_no_leap


def take_the_polygon_with_largest_area_if_mutlipolygon(intersection_area):
    polygon = max(intersection_area,
                  key=lambda a: a.area) if intersection_area.geom_type == 'MultiPolygon' else intersection_area
    return polygon


def next_month(month):
    if month == 11:
        new_month = 0
    else:
        new_month = month+1
    return new_month


def previous_month(month):
    if month == 0:
        new_month = 11
    else:
        new_month = month-1
    return new_month


def fwd_days_to_next_month(start_day=0, window_size=20, month=0, year=0):
    old_month = month
    old_year = year
    days_in_this_month = leap_day(year)[month_names[month]]

    year_month_day_index = []
    year_month_day_name = []

    if ((start_day+window_size) >= days_in_this_month):
        new_month = next_month(month)
        new_year = year + \
            1 if ((old_month == 11) and ((new_month == 0)
                                         or (new_month == 1))) else old_year

    for i in range(window_size):
        day = start_day+i
        if day >= days_in_this_month:
            day = day-days_in_this_month
            month = new_month
            year = new_year
        else:
            month = old_month
            year = old_year
        year_month_day_index.append((year, month, day))
        year_month_day_name.append(
            ('%d-%s-%d' % (year+ZEROTH_YEAR, month_names[month], day+1)))
    return year_month_day_index, year_month_day_name


# logging_object.write('DJF, MAM, JJA, SON')
# DJF = [10, 11, 0, 1]
# MAM = [1, 2, 3, 4]
# JJA = [4, 5, 6, 7]
# SON = [7, 8, 9, 10]

def return_window_dates(month_index=[10,11,0,1], orig_year=0, window_size=20, starting_day_for_previous_month=15):
    indices_each_year = []
    date_each_year = []
    for mo in month_index:
        year = orig_year-1 if ((mo == 10) or (mo == 11)) else orig_year
        starting_day = starting_day_for_previous_month if mo == month_index[0] else 0
        for D in range(starting_day, leap_day(year)[month_names[mo]]):
            indices_each_year.append(fwd_days_to_next_month(
                start_day=D, window_size=window_size, month=mo, year=year)[0])
            date_each_year.append(fwd_days_to_next_month(
                start_day=D, window_size=window_size, month=mo, year=year)[1])
    return indices_each_year, date_each_year


# logging_object.write('Defined more functions')


def bwd_days_to_previous_month(start_day=0, window_size=20, month=0, year=0):
    old_month = month
    old_year = year
    days_in_this_month = leap_day(year)[month_names[month]]

    year_month_day_index = []
    year_month_day_name = []

    if ((start_day-window_size) < 0):
        new_month = previous_month(month)
        new_year = year - \
            1 if ((old_month == 0) and ((new_month == 11)
                                        or (new_month == 10))) else old_year
        days_in_this_month = leap_day(year)[month_names[new_month]]

    for i in range(window_size):
        day = start_day-i
        if day < 0:
            day = days_in_this_month+day
            month = new_month
            year = new_year
        else:
            month = old_month
            year = old_year
        year_month_day_index.append((year, month, day))
        year_month_day_name.append(
            ('%d-%s-%d' % (year+ZEROTH_YEAR, month_names[month], day+1)))
    return year_month_day_index, year_month_day_name


def plus_minus_10_days(start_day=0, window_size=20, month=0, year=0):
    plus_days = fwd_days_to_next_month(start_day, window_size, month, year)[1]
    minus_days = bwd_days_to_previous_month(
        start_day, window_size, month, year)[1][1:][::-1]

    plus_days_tuple = fwd_days_to_next_month(
        start_day, window_size, month, year)[0]
    minus_days_tuple = bwd_days_to_previous_month(
        start_day, window_size, month, year)[0][1:][::-1]

    all_days = minus_days + plus_days
    all_days_tuple = minus_days_tuple + plus_days_tuple

    return all_days, all_days_tuple



def contour_coordinates(field, critical=55, month=0, year=0, day=0):
#     field = A_N_daily
    def add_cyclic(YY, x=None):
        x = lon
        YY = np.append(YY, YY[..., :int(len(x)/2)], axis=-1)
        return YY

    def latavg(x, year=0, month=0, day=0):
        x0 = x[year, month, day, ...]
        x1 = x0[::2, :]
        x2 = x0[1::2, :]
        x = (x1+x2)*0.5
        return x[:-1]

    lon_extended = np.append(lon, lon[:int(len(lon)/2)]+360)

    if isinstance(critical, (list, tuple, np.ndarray)):
        critical = np.copy(critical)[month, ...]

    field = latavg(np.copy(field), year, month, day)

    month_day_year = ('%2d/%2d/%4d' % (month+1, day+1, year+ZEROTH_YEAR))

    xopen = []
    xclose = []
    centroids = []
    centroids_shifted = []

    def get_contours(array=add_cyclic(field-critical), threshold=0):
        contour_points = []
        contours = measure.find_contours(array, threshold)
        x = lon_extended
        y = lat_mean
        for k in range(0, len(contours)):
            contour = contours[k]
            contour_x = x[0] + contour[:, 1] * (x[-1] - x[0]) / (len(x)-1)
            contour_y = y[0] + contour[:, 0] * (y[-1] - y[0]) / (len(y)-1)
            contour_points.append((contour_x, contour_y))
        return contour_points

    contours = get_contours(add_cyclic(field-critical), 0)

    def close_some_necessary_open_contours(x, y):
        xmid = np.mean([x[0], x[-1]])
        ymid = np.mean([y[0], y[-1]])
        xnew = np.append(x, xmid)
        ynew = np.append(y, ymid)
        xclosed = np.append(xnew, x[0])
        yclosed = np.append(ynew, y[0])
        return (xclosed, yclosed)

    def check_closed_contours(x, y, ii):
        if (y[0] != y[-1]):
            lab = 'x open '+str(ii)
            xopen.append((x, y))
        else:
            lab = 'x close '+str(ii)
            if (x[0] != x[-1]):  # these are the open contours at lower and upper latitudinal bounds
                x, y = close_some_necessary_open_contours(x, y)
            centroids.append(centroid(x, y))
            if centroid(x, y)[0] >= 180:
                centroids_shifted.append(
                    (centroid(x, y)[0]-360, centroid(x, y)[1]))
            else:
                centroids_shifted.append(centroid(x, y))
            xclose.append((x, y))
        return lab

    def compare_coord(coords1, coords_array):
        for index, values in enumerate(coords_array):
            if ((np.abs(coords1[0]-values[0]) < 1e-6) & (np.abs(coords1[1]-values[1]) < 1e-6)):
                return False
        return True

    def unwrap(values):
        val = np.copy(values)
        val[0][val[0] > 180] = val[0][val[0] > 180]-360
        return val

    def avg_lat_lon_extent(x, y):
        xcoords = np.copy(x)
        ycoords = np.copy(y)
        x_max = np.max(np.abs(xcoords-np.mean(xcoords)))
        y_max = np.max(np.abs(ycoords-np.mean(ycoords)))
        return x_max, y_max

    lab                         = [check_closed_contours(seg[0], seg[1], ii) for ii, seg in enumerate(contours)]                     
    centroids_remove_duplicates = [values for num, values in enumerate(centroids_shifted) if compare_coord(values, centroids_shifted[:num])]                    
    xclose_remove_duplicates    = [xclose[num] for num, values in enumerate(centroids_shifted) if compare_coord(values, centroids_shifted[:num])]
    xclose_unwrap_contours      = [unwrap(values) for num, values in enumerate(xclose_remove_duplicates)]
    areas                       = [area(seg[0], seg[1]) for ii,     seg in enumerate(xclose_remove_duplicates)]
    xyrange                     = [avg_lat_lon_extent(seg[0], seg[1]) for ii,     seg in enumerate(xclose_remove_duplicates)]

    def convert_to_shapely_objects(contours=xclose_remove_duplicates):

        def return_tuple(contours, k, i):
            x = contours[k][0][i]
            y = contours[k][1][i]
            return (x, y)

        def transform(coords):
            R_sq = (6371e3)**2
            x = np.deg2rad(coords[0])*R_sq
            y = np.sin(np.deg2rad(coords[1]))
            return (x, y)
        shapely_contours           = shp.MultiPolygon([shp.Polygon([return_tuple(contours, k, i) for i in range(len(contours[k][0]))]) for k in range(len(contours))])
        shapely_area               = [(shp.Polygon([transform(return_tuple(contours, k, i)) for i in range(len(contours[k][0]))])).area for k in range(len(contours))]
        # This area is using transformed coordinates, so values are in m^2
        shapely_centroid_unbounded = ([(shp.Polygon([return_tuple(contours, k, i) for i in range(len(contours[k][0]))])).centroid for k in range(len(contours))])
        shapely_centroid_bounded   = ([unwrap(values.coords)[0] for values in shapely_centroid_unbounded])
        shapely_dict               = {'contours': shapely_contours, \
                                      'area': shapely_area, \
                                      'shapely_centroid_unbounded': shp.MultiPoint(shapely_centroid_unbounded),
                                      'shapely_centroid_bounded': shp.MultiPoint(shapely_centroid_bounded)}
        return shapely_dict

    shapely_dict = convert_to_shapely_objects(xclose_remove_duplicates)

    # 1. Contours --> Has all the contours. There are some redundant ones which needs to be removed.
    # 2. xclose   --> Has only the closed contours. The x-open contours, i.e, contours that are open along the latitude axes are discarded. But has many redundants
    # 3. xclose_remove_duplicates --> x-close contours but with the redundant ones removed.
    # 4. xclose_unwrap_contours --> The same x-close contours as (3), but the points where x>180, the x coordinates are shifted back to the range -180 to 180.
    # 5. centroids --> Centroids of all the contours. Has many redundants
    # 6. centroids_shifted --> Same as (5) but the x coordinate range is kept between -180 to 180 but it still has redundants
    # 7. centroids_remove_duplicates --> Sames as (6) but with duplicate centroids removed.
    # 8. area --> Gives the area of enclosed by the contours.
    # 9. xy range

    contour_data = {'lon': lon, 'lat': lat, \
                    'date': month_day_year, \
                    'contours': xclose_unwrap_contours,  \
                    'centroids': centroids_remove_duplicates,
                    'raw_area': areas, \
                    'xyrange': xyrange, \
                    'contours_unshifted': xclose_remove_duplicates, \
                    'shapely_dict': shapely_dict}

    # The Mth coordinate of the Nth contour can be accessed by xclose[N][M][x,y]
    # The centroid of the Nth contour can be accessed by centroid[N][x,y]

    return contour_data


def plot_contour(YY, ms='o', markersize=20, mfc='red', lw=1):
    for n in range(len(YY)):
        xcoord = ((YY)[n])[0]
        ycoord = ((YY)[n])[1]
        py.plot(xcoord, ycoord, ms, markersize=markersize,
                label=str(n), mfc=mfc, lw=lw)


def plot(combine, alpha=0.5, labels='', color='r'):
    norm = py.cm.colors.Normalize(vmin=0, vmax=len(combine))
    for poly0 in (combine):
        py.plot(*poly0.exterior.xy, color='k', ls='-', alpha=0.5)
        py.fill(*poly0.exterior.xy, color=color, ls='-', alpha=alpha)
        py.text(*poly0.centroid.coords[0], s=str(labels), fontsize=20)


def block_days(BLOCKS, co, N):
    alphas = np.linspace(0.5, 0.1, len(BLOCKS))
    for ii, day in enumerate(BLOCKS):
        plot(day, alpha=alphas[ii], color=co, labels=N)
    py.xlim(-180, 250)
    py.ylim(20, 70)


def intersection_is_significant(intersection, contour1, contour2, threshold=0.2):

    if isinstance(intersection, list):
        intersected_area = np.sum([inter.area for inter in intersection])
    else:
        intersected_area = intersection.area
    fraction = intersected_area/np.mean([contour1.area, contour2.area])
    if (fraction > threshold):
        return True
    else:
        return False


def search_across_contours(contours1, contours2, threshold=0.2):
    comb0 = []
    comb1 = []
    for contour1 in contours1:
        for contour2 in contours2:
            intersection_exists = take_the_polygon_with_largest_area_if_mutlipolygon(
                contour1.intersection(contour2))
            if ((intersection_exists) and (intersection_is_significant(intersection_exists, contour1, contour2, threshold))):
                # Collects all the contours of day 1 which has intersections with day 2.
                comb0.append(contour1)
                # Collects all the contours of day 2 which has intersections with day 1.
                comb1.append(contour2)
    return comb0, comb1


def minus(superset, subset):
    result = [x for x in superset if x not in subset]
    return result


def create_blocking_matrix(window_size_days=6, threshold=0.2, field=None, critical=55,
                           year=0, start_day=0, month_index=[10, 11, 0, 1], starting_day_for_previous_month=15):

    ## field = A_N_daily
    M = [[np.nan] * (window_size_days) for i in range(window_size_days)]

    date_tuple = return_window_dates(month_index=month_index, orig_year=year, window_size=window_size_days,\
                                     starting_day_for_previous_month=starting_day_for_previous_month)[0][start_day]  # Creates a whole list of
    date_names = return_window_dates(month_index=month_index, orig_year=year, window_size=window_size_days,
                                     starting_day_for_previous_month=starting_day_for_previous_month)[1][start_day]

    # Fill the first row of the matrix
    start_time = ti.time()
    for i in range(window_size_days):
        if ((start_day > 0)):
            date_indices_minus1 = return_window_dates(
                month_index=month_index, orig_year=year, window_size=window_size_days)[0][start_day-1]
            if (i == 0):
                #  day  =  start_day+i
                y, m, d = date_tuple[i]
                all_contours_on_start_day_dict = (contour_coordinates(field=field, critical=critical, month=m, year=y, day=d)['shapely_dict'])
                all_contours_on_start_day      = [(all_contours_on_start_day_dict['contours'].geoms[N]) for N in range(len(all_contours_on_start_day_dict['contours']))]

                y, m, d = date_indices_minus1[0]
                all_contours_previous_day_dict = (contour_coordinates(field=field, critical=critical, month=m, year=y, day=d)['shapely_dict'])
                all_contours_previous_day      = [(all_contours_previous_day_dict['contours'].geoms[N]) for N in range(len(all_contours_previous_day_dict['contours']))]

                intersected_contours_between_the_two_days, discard = search_across_contours(contours1=all_contours_on_start_day, contours2=all_contours_previous_day, threshold=threshold)
                contours_needed_for_start_day = minus(all_contours_on_start_day, intersected_contours_between_the_two_days)

                M[0][i] = contours_needed_for_start_day

            else:
                # day = start_day+i
                y, m, d = date_tuple[i]
                shapely_dict = (contour_coordinates(field=field, critical=critical, month=m, year=y, day=d)['shapely_dict'])
                contour_day  = [(shapely_dict['contours'].geoms[N]) for N in range(len(shapely_dict['contours']))]
                M[0][i]      = contour_day
                #print (i, end=',')

        if start_day == 0:
            # day  = start_day+i
            y, m, d = date_tuple[i]
            shapely_dict = (contour_coordinates(field=field, critical=critical, month=m, year=y, day=d)['shapely_dict'])
            contour_day  = [(shapely_dict['contours'].geoms[N]) for N in range(len(shapely_dict['contours']))]
            M[0][i] = contour_day
            #print (i, end=',')

    end_time = ti.time()
    time1 = (end_time - start_time)

    # Create a matrix of blocks starting day 0

    def adjust_i1(j):
        if j == 0:
            return i+1
        else:
            return i+2

    def adjust_i2(i):
        if i == 0:
            return i
        else:
            return i+1

    def Ms(i, j):
        print('M[%d][%d]' % (adjust_i1(j), j), 'M[%d][%d]' % (i+1, j+1), ' = M[%2d][%2d]' % (adjust_i1(j)-1, j),
              'M[%2d][%2d]' % (adjust_i2(i), j+1))  # Just for sanity check of the matrix sequence elements

    start_time = ti.time()
    for i in range(0, window_size_days-1):
        for j in range(0, window_size_days-1-i):
            #             Ms(i,j)
            M[adjust_i1(j)][j], M[i+1][j+1] = search_across_contours(contours1=M[adjust_i1(j)-1][j], contours2=M[adjust_i2(i)][j+1], threshold=threshold)
#         print ('--------------')
    end_time = ti.time()
    time2 = (end_time - start_time)

    return time1, time2, M, date_names, date_tuple


# logging_object.write('Defined some more functions')


def plot_nth_day(start_day=0, field=None, critical=55, year=0,
                 fill=True, plot_boolean=True, window_size_days=20, threshold=0.4, month_index=[10,11,0,1],
                 starting_day_for_previous_month=15):

    ## field = A_N_daily
    time0 = []
    time1 = []

    block_matrix_and_time = create_blocking_matrix(window_size_days=window_size_days, threshold=threshold,
                                                   start_day=start_day, year=year, field=field,
                                                   critical=critical, month_index=month_index,
                                                   starting_day_for_previous_month=starting_day_for_previous_month)
    time0.append(block_matrix_and_time[0])
    time1.append(block_matrix_and_time[1])
    M          = block_matrix_and_time[2]
    date_names = block_matrix_and_time[3]  # Date on string format
    # Date in tuple format (year, month, day) ###
    date_tuple = block_matrix_and_time[4]

    # It only has length of the window size since start day is already taken
    # This returns all the blocks that emerged on this particular date.

    def plot(combine, alpha=0.5, labels='', color='r'):
        norm = py.cm.colors.Normalize(vmin=0, vmax=len(combine))
        for poly0 in (combine):
            py.plot(*poly0.exterior.xy, color='k', ls='-', alpha=0.5)
            if fill:
                py.fill(*poly0.exterior.xy, color=color, ls='-', alpha=alpha)
            py.text(*poly0.centroid.coords[0], s=str(labels), fontsize=20)

    def block_days(BLOCKS, co, N):
        alphas = np.linspace(0.5, 0.1, len(BLOCKS))
        for ii, day in enumerate(BLOCKS):
            plot(day, alpha=alphas[ii], color=co, labels=N)
        py.xlim(-180, 250)
        py.ylim(20, 70)

    block_day_contours = [[] for i in range(window_size_days)]

    # N day blocked_contour
    def sumj(j, window=window_size_days):
        i = window-1 if j == 0 else window-j
        return i

    def fill_blocks(N):
        block_day_contours[N] = [
            minus(M[sumj(j, N+2)-1][j], M[sumj(j, N+2)][j]) for j in range(N+1)]

    block_day_contours[-1] = [M[sumj(j)][j] for j in range(window_size_days)]
    for N in range(0, window_size_days-1):
        fill_blocks(N)

    if plot_boolean:
        py.figure(figsize=(30, 5))
        co = py.cm.hsv(np.linspace(0, 1, window_size_days))
        for N in range(1, window_size_days+1):
            dayN = block_day_contours[N-1]
            block_days(dayN, co[N-1], N)
            #print (N, end=',')

    return block_day_contours, date_names, date_tuple


def not_empty(list1):
    qual = 0 if len(list1) == 0 else 1
    return qual


def intersect_multiple_polygons(circles):
    intersections = []
    idx = index.Index()

    for pos, circle in enumerate(circles):
        idx.insert(pos, circle.bounds)

    for circle in circles:
        merged_circles = ops.unary_union([circles[pos] for pos in idx.intersection(circle.bounds) if circles[pos] != circle])
        intersections.append(circle.intersection(merged_circles))
    intersection = ops.unary_union(intersections)
    return intersection


def xy_coords_separately(C1):
    xcoords = [c[0] for c in C1]
    ycoords = [c[1] for c in C1]
    coords = (xcoords, ycoords)
    return coords


def contours_containing_point(centroidi, polygonsi):
    stationary_days = []
    for ii, poly in enumerate(polygonsi):
        point = shp.Point(centroidi[0])
        if point.within(poly):
            stationary_days.append(ii)
    return stationary_days


def polar_area(geom):
    geom_area = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init='EPSG:4326'),
            pyproj.Proj(
                proj='aea',
                lat_1=geom.bounds[1],
                lat_2=geom.bounds[3])),
        geom)
    return (geom_area.area)/10**6  # return in km ^2
    # This area is using transformed coordinates, so values are in m^2


def latavg(x, year=0, month=0, day=0):
    x0 = x[year, month, day, ...]
    x1 = x0[::2, :]
    x2 = x0[1::2, :]
    x = (x1+x2)*0.5
    return x[:-1]


def time_evolution_of_blocks_at_the_centroid(field=None, centroid=(0, 0), lon=None, lat=None, days_around_peak_blocking=[(1, 0, 0)],
                                             evaluating_critical_field=False, do_lat_avg=True):

    ## field = A_N_daily
    #### Right now lat_mean is only taking for NH ####
    field_at_loc = []

    lon_at, lat_at = centroid
    if lon_at > 180:
        lon_at = lon_at-360

    diff_la = np.abs(lat-lat_at)
    diff_lo = np.abs(lon-lon_at)

    ind_la = np.squeeze(np.where(diff_la == np.min(diff_la)))
    ind_lo = np.squeeze(np.where(diff_lo == np.min(diff_lo)))

    for dates in days_around_peak_blocking:
        y, m, d = dates
#        logging_object.write('time evolution : (%d/%d/%d)'%(y, m, d))

        if evaluating_critical_field:
            field1 = np.copy(field)[m, ...]

        else:
            if (do_lat_avg):
                field1 = latavg(np.copy(field), y, m, d)
            else:
                field1 = np.copy(field)[y, m, d, ...]

        field_at_loc.append(field1[ind_la, ind_lo])

    return field_at_loc


def latavg_without_ymd(x):
    x0 = x
    x1 = x0[..., ::2, :]
    x2 = x0[..., 1::2, :]
    x = (x1+x2)*0.5
    return x[..., :-1, :]


import itertools
def same_or_near_day(X1,X2, block_length):
    from datetime import date
    
    X1=X1[0]; X2=X2[0]
    def year_of(Y):
        return int(Y.split('-')[0])
    def month_of(Y):
        return month_index[Y.split('-')[1]]
    def day_of(Y):
        return int(Y.split('-')[2])
        
    y1_date       =  date(year_of(X1), month_of(X1), day_of(X1))
    y2_date       =  date(year_of(X2), month_of(X2), day_of(X2))
    delta_days    =  np.abs((y2_date - y1_date).days)
    if delta_days<block_length:
        return True
    else:
        return False
       
def not_very_distant(X1,X2):    
    p0, p1 = X1[0][0], X2[0][0]
    latitude_distance  = np.abs(p1[1]-p0[1])
    longitude_distance = np.abs(p1[0]-p0[0])
    RETURN=True if ((latitude_distance<4) and (longitude_distance<20)) else False
    return RETURN
                 
def compare(keys_pair, block_year):
    key1 = keys_pair[0]
    key2 = keys_pair[1]
    max_block_length = max(block_year[key1]['avg_no_of_blocked_days'][0], block_year[key2]['avg_no_of_blocked_days'][0])
    time_comparison  = same_or_near_day(block_year[key1]['peak_block_date'],   block_year[key2]['peak_block_date'], max_block_length)
    space_comparison = not_very_distant(block_year[key1]['envelope_centroid'], block_year[key2]['envelope_centroid'])
    return (time_comparison and space_comparison)    

def minimum_key(keys_pair, block_year):
    key1 = keys_pair[0]
    key2 = keys_pair[1]
    if (block_year[key1]['avg_no_of_blocked_days'][0] < block_year[key2]['avg_no_of_blocked_days'][0]):
        deletion_key = key1
    else:
        deletion_key = key2
    return deletion_key
    

def remove_redundant_blocks(BLOCKS_IN_A_YEAR):
    deletion_keys=[]
    for keys_pair in itertools.combinations(BLOCKS_IN_A_YEAR.keys(), 2):
        if compare(keys_pair, BLOCKS_IN_A_YEAR):
            deletion_key=minimum_key(keys_pair, BLOCKS_IN_A_YEAR)
            deletion_keys.append(deletion_key)
    for deletion_key in deletion_keys:
        BLOCKS_IN_A_YEAR.pop(deletion_key, None)
        
    return BLOCKS_IN_A_YEAR


def save_blocking_info(start_day, field=None, critical=55, fill=True, plot_boolean=True,
                       window_size_days=12, threshold=0.6, year=2, starting_day_for_previous_month=20, month_index=[10,11,0,1],
                       lon=None, lat=None, evolve_window=12, input_dict={}, save_separate=False, DEST=None):

    ### field = A_N_daily
    start_N = ti.time()
    logging_object.write('--------------------------------------')
    logging_object.write('N= %d' % (start_day))
    logging_object.write('--------------------------------------')

    category_of_contours_starting_day_A, date_names, date_tuple = plot_nth_day(start_day=start_day, field=field,
                                                                               critical=critical, fill=fill, plot_boolean=plot_boolean,
                                                                               window_size_days=window_size_days, threshold=threshold,
                                                                               starting_day_for_previous_month=starting_day_for_previous_month,
                                                                               year=year, month_index=month_index)

    logging_object.write('Contours extracted by creating block matrix')

    # Denotes the starting day when blocks have emerged.
    y_start, m_start, d_start = date_tuple[0]

    logging_object.write('<------ Start date is %s ------->' % (date_names[0]))

    C1 = []
    cxy1 = []
    I1 = []
    Ic1 = []
    max_days1 = []
    stationary_daysI1 = []
    avg_daysI1 = []
    total_area1 = []
    avg_area1 = []
    areas = []
    max_block_date = []
    peak_block_date = []
    stationary_dates_I1 = []
    peak_block_day = []
    blocking_A = []
    blocking_A1 = []
    blocking_A2 = []
    blocking_U = []
    blocking_U_URB = []
    blocking_dates = []
    blocking_time = []
    blocking_Ac = []
    blocking_U_U0 = []
    blocking_F = []
    blocking_F1 = []
    blocking_F2 = []
    blocking_F3 = []
    latitude_bands = []

    logging_object.write('Created all empty arrays')

    for i in range(len(category_of_contours_starting_day_A)):

        logging_object.write('Sweeping through day %d starting %s' % (i, date_names[0]))

        co = py.cm.hsv(np.linspace(0, 1, len(category_of_contours_starting_day_A)))
        # starting day 0 but last for i+1 days
        Dayi = category_of_contours_starting_day_A[i]

        if (not_empty(Dayi[0]) and (i > 1)):
            min_number_of_polygon_set = np.min([len(Dayi[i]) for i in range(len(Dayi))])
            polygonssi = [[poly[k] for poly in Dayi] for k in range(min_number_of_polygon_set)]

            for polygonsi in polygonssi:

                # date_names[i] --> denotes the end date of a certain polygon (Note that 0<=i<20)

                intersectioni = take_the_polygon_with_largest_area_if_mutlipolygon(intersect_multiple_polygons(polygonsi))
                centroid_intersecti = intersectioni.centroid.coords
                blocking_days = contours_containing_point(centroid_intersecti, polygonsi)

                x, y = centroid_intersecti[0]
                logging_object.write('Centroid lon,lat = (%1.2f, %1.2f)' % (x, y))

                if y < 30:
                    lat_band = '20-30'
                elif ((y >= 30) and (y < 45)):
                    lat_band = '30-45'
                elif ((y >= 45) and (y < 60)):
                    lat_band = '45-60'
                elif ((y >= 60) and (y < 75)):
                    lat_band = '60-75'
                else:
                    lat_band = '75-90'

                if len(blocking_days) >= 1:

                    logging_object.write('Now on day %s for blocks> 1 days' %(date_names[i]))

                    latitude_bands     .append(lat_band)
                    # All the contours of a given length
                    C1                 .append(polygonsi)

                    # All the centroids saved as [(xi,yi)]
                    c1 = [poly.centroid.coords[0] for poly in polygonsi]

                    # All the centroids saved as ([xi],[yi])
                    cxy1               .append(xy_coords_separately(c1))

                    I1                 .append(intersectioni)
                    Ic1                .append(centroid_intersecti)

                    max_days1          .append(len(polygonsi))
                    max_block_date     .append(date_names[i])

                    stationary_daysI1  .append(blocking_days)
                    stationary_dates_I1.append([date_names[int(n)] for n in blocking_days])

                    peak_block_day     .append(blocking_days[int(len(blocking_days)/2)])
                    peak_block_date    .append(date_names[int(blocking_days[int(len(blocking_days)/2)])])

                    avg_daysI1         .append(len(blocking_days))
                    # Gives a measure of the days the centroid of the union has high wave activity

                    total_area1        .append(polar_area(intersectioni))
                    avg_area1          .append(np.mean(np.array([polar_area(poly) for poly in polygonsi])))
                    areas              .append([polar_area(poly) for poly in polygonsi])

                    y_max,  m_max,   d_max = date_tuple[i]
                    y_peak, m_peak,  d_peak = date_tuple[int(blocking_days[int(len(blocking_days)/2)])]

                    plus_minus = plus_minus_10_days(start_day=d_peak, window_size=evolve_window, month=m_peak, year=y_peak)
                    days_around_peak_blocking  = plus_minus[1]
                    dates_around_peak_blocking = plus_minus[0]

                    logging_object.write('Appended some data')

                    logging_object.write('--- LWA ----')
                    start = ti.time()
                    field_A = time_evolution_of_blocks_at_the_centroid(field=A_N_daily, centroid=centroid_intersecti[0],
                                                                       lon=lon, lat=lat_mean,  days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of LWA : %1.2f' % (end-start))

                    logging_object.write('---U ----')
                    start = ti.time()
                    field_U = time_evolution_of_blocks_at_the_centroid(field=U_N_daily,  centroid=centroid_intersecti[0],
                                                                       lon=lon, lat=lat_mean,  days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of U : %1.2f' % (end-start))

                    logging_object.write('---U/URB ----')
                    start = ti.time()
                    field_U_URB = time_evolution_of_blocks_at_the_centroid(field=Uratio_N_daily, centroid=centroid_intersecti[0],
                                                                           lon=lon, lat=lat_mean, days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of U/URB: %1.2f' % (end-start))


                    logging_object.write('---F ----')
                    start = ti.time()
                    field_F = time_evolution_of_blocks_at_the_centroid(field=F_N_daily,  centroid=centroid_intersecti[0], lon=lon, lat=lat_mean,
                                                                       days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of F: %1.2f' % (end-start))

                    logging_object.write('---F2 ----')
                    start = ti.time()
                    field_F1 = time_evolution_of_blocks_at_the_centroid(field=F2_N_daily,  centroid=centroid_intersecti[0], lon=lon, lat=lat_mean,
                                                                        days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of F2: %1.2f' % (end-start))

                    logging_object.write('---F1 ----')
                    start = ti.time()
                    field_F2 = time_evolution_of_blocks_at_the_centroid(field=F1_N_daily,  centroid=centroid_intersecti[0], lon=lon, lat=lat_mean,
                                                                        days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of F1: %1.2f' % (end-start))

                    logging_object.write('---F3 ----')
                    start = ti.time()
                    field_F3 = time_evolution_of_blocks_at_the_centroid(field=F3_N_daily,  centroid=centroid_intersecti[0], lon=lon, lat=lat_mean,
                                                                        days_around_peak_blocking=days_around_peak_blocking)
                    end = ti.time()
                    logging_object.write('Time evolution of F3: %1.2f' % (end-start))

                    logging_object.write('---Critical Ac ----')
                    # If evaluating critical field, do_lat_avg doesnt matter
                    start = ti.time()
                    if isinstance(critical, (list, tuple, np.ndarray)):

                        field_critical_Ac = time_evolution_of_blocks_at_the_centroid(field=critical,  centroid=centroid_intersecti[0], lon=lon, lat=lat_mean,
                                                                                     days_around_peak_blocking=days_around_peak_blocking, evaluating_critical_field=True)
                    else:
                        field_critical_Ac = critical * (np.ones(len(days_around_peak_blocking)))
                    end = ti.time()
                    logging_object.write('Time evolution of critical Ac: %1.2f' %(end-start))
                    
                    
                    # Need to write a few more lines of code to save F, F1, F2, F3

                    blocking_A.append(field_A)
                    
                    blocking_F.append(field_F)
                    blocking_F1.append(field_F1)
                    blocking_F2.append(field_F2)
                    blocking_F3.append(field_F3)
                    blocking_U.append(field_U)
                    blocking_U_URB.append(field_U_URB)
                    blocking_Ac.append(field_critical_Ac)
                    blocking_time.append(days_around_peak_blocking)
                    blocking_dates.append(dates_around_peak_blocking)
            
                    logging_object.write('Appended time evolution of blocking events')

    # envelope centroid is the block point chosen

    if not(not(avg_daysI1)):

        blocking_info = {'contours': C1, 'centroids': cxy1,    'min_2day_envelope': I1,   'envelope_centroid': Ic1, \
                         'max_blocked_days': max_days1,        'blocked_days_of_centroid': stationary_daysI1, \
                         'avg_no_of_blocked_days': avg_daysI1, 'total_area_of_blocking_event': total_area1, \
                         'peak_block_day': peak_block_day,     'avg_area_of_blocking_event': avg_area1, 'areas': areas, \
                         'peak_block_date': peak_block_date,   'blocked_dates_of_centroid': stationary_dates_I1, \
                         'max_block_date': max_block_date,     'blocking_A': blocking_A, 'blocking_time': blocking_time, \
                         'blocking_dates': blocking_dates,     'start_block_day': date_tuple[0], 'start_block_date': date_names[0], \
                         'blocking_F': blocking_F,             'blocking_F1': blocking_F1,       'blocking_F2': blocking_F2, \
                         'blocking_F3': blocking_F3,           'blocking_U': blocking_U,\
                         'blocking_Ac': blocking_Ac,           'latitude_bands': latitude_bands}
        
#                          'blocking_F': blocking_F, 'blocking_F1': blocking_F1, 'blocking_F2': blocking_F2, 'blocking_F3': blocking_F3, 'blocking_U': blocking_U, \

                                                                                             
        save_data = {'start_block_date='+date_names[0]: blocking_info}

        if save_separate:
            DEST = DEST+'/%d/' % (THIS_YEAR+ZEROTH_YEAR)
            h5saveload.make_sure_path_exists(DEST)
            ### '/data/pragallva/Spring_quarter_2020/post_processing/blocking_info/%s/Ac=55/%d/' % ('NH',THIS_YEAR+ZEROTH_YEAR)
#             h5saveload.save_dict_to_hdf5(blocking_info, DEST+date_names[0]+'.hdf5',)
            hkl.dump(blocking_info, DEST+date_names[0]+'.hkl',)

        input_dict.update(save_data)

        logging_object.write('Updated dictionary')

        if plot_boolean:
            py.title(blocking_info['start_block_date'], fontsize=30)

    else:
        logging_object.write('No blocks detected')

    end_N = ti.time()
    logging_object.write('FINISHED ALL SWEEP FOR N=%d' % (start_day))
    logging_object.write('Total time taken %1.2f' % (end_N-start_N))
    logging_object.write('--------------------------------------')

    return input_dict, date_names[0]



###################### Defined some more functions to remove remaining duplicate events ##################

def return_sorted_blocks(year, blocking_data):
    start=ti.time()
    year   = '%d.hkl'%(year)
    blocks =  hkl.load(blocking_data+year)
    end=ti.time()
    block_keys     = list(blocks.keys())

    start_date_object = []
    for i in range(len(block_keys)):
        y,m,d=blocks[block_keys[i]]['start_block_day']
        y =y+1979; m=m+1; d=d+1
        datetime_object = datetime(y,m,d)
        blocks[block_keys[i]]['start_block_day']=datetime_object

    new_blocks={}
    listofTuples=sorted(blocks.items(), key = lambda x: x[1]['start_block_day'])
    for elem in listofTuples :
        new_blocks.update({elem[0]: elem[1]})

    end=ti.time()
    print (year, 'time taken = %1.2f'%(end-start) ) #, end=', ')
    return new_blocks


def same_or_near_day_all_blocks(X1,X2, block_length):
    from datetime import date        
#     X1=X1[0]; X2=X2[0]
    def year_of(Y):
        return int(Y.split('-')[0])
    def month_of(Y):
        return month_index[Y.split('-')[1]]
    def day_of(Y):
        return int(Y.split('-')[2])
        
    y1_date       =  date(year_of(X1), month_of(X1), day_of(X1))
    y2_date       =  date(year_of(X2), month_of(X2), day_of(X2))
    delta_days    =  np.abs((y2_date - y1_date).days)
    if delta_days < block_length:
        return True
    else:
        return False
    
def not_very_distant_all_blocks(X1,X2):    
#     p0, p1 = X1[0][0], X2[0][0]
    p0, p1 = X1[0], X2[0]
    latitude_distance  = np.abs(p1[1]-p0[1])
    longitude_distance = np.abs(p1[0]-p0[0])
    RETURN=True if ((latitude_distance<4) and (longitude_distance<10)) else False
    return RETURN

def minimum_key_all_blocks(keys_pair, block_year, j,k):
    key1 = keys_pair[0]
    key2 = keys_pair[1]
    if (block_year[key1]['avg_no_of_blocked_days'][j] < block_year[key2]['avg_no_of_blocked_days'][k]):
        deletion_key,list_index = (key1,j)
    else:
        deletion_key,list_index = (key2,k)
    return deletion_key,list_index

                 
def compare_and_delete(keys_pair, block_year,DELETE_INDEX={}, year=1980):
        
    key1 = keys_pair[0]
    key2 = keys_pair[1]
         
    for j in range(len(block_year[key1]['avg_no_of_blocked_days'])):
        for k in range(len(block_year[key2]['avg_no_of_blocked_days'])):
            max_block_length = max(block_year[key1]['avg_no_of_blocked_days'][j], block_year[key2]['avg_no_of_blocked_days'][k])
            time_comparison  = same_or_near_day_all_blocks(block_year[key1]['peak_block_date'][j],   block_year[key2]['peak_block_date'][k], max_block_length)
            space_comparison = not_very_distant_all_blocks(block_year[key1]['envelope_centroid'][j], block_year[key2]['envelope_centroid'][k])

            ### delete it here itself
            if (time_comparison and space_comparison):
                deletion_key,list_index = minimum_key_all_blocks(keys_pair, block_year, j,k)
                all_keys=list(block_year[deletion_key].keys())
                all_keys.remove('start_block_date')
                all_keys.remove('start_block_day')
                
                if deletion_key in block_year:
                    for every_key in all_keys:                    
                        if every_key in block_year[deletion_key]: #.has_key(every_key):                           
                             DELETE_INDEX.update({deletion_key:[year,list_index]})
                                
    return DELETE_INDEX                            
                                



################## All the input parameters  ##################


if __name__ == "__main__":
    
        current_file_name = sys.argv[0].split('/')[-1].split('.py')[0]
        today = date.today()
        cur_date = today.strftime("%d.%m.%Y")
        cur_file_date  = '%s-%s'%(current_file_name, cur_date)

    
        begin_year  = 1980 #1980  ## int(sys.argv[1]) ##
        finish_year = 2022 #2022  ### int(sys.argv[2]) ##
        
        month_names  = {0: 'Jan',   1: 'Feb',    2: 'Mar',    3: 'Apr',    4: 'May',   5: 'Jun',
                        6: 'Jul',    7: 'Aug',    8: 'Sep',    9: 'Oct',    10: 'Nov',    11: 'Dec'}
        days_no_leap = {'Jan': 31,  'Feb': 28,   'Mar': 31,   'Apr': 30,   'May': 31,  'Jun': 30,
                        'Jul': 31,   'Aug': 31,   'Sep': 30,   'Oct': 31,    'Nov': 30,    'Dec': 31}
        days_leap    = {'Jan': 31,  'Feb': 29,   'Mar': 31,   'Apr': 30,   'May': 31,  'Jun': 30,
                        'Jul': 31,   'Aug': 31,   'Sep': 30,   'Oct': 31,    'Nov': 30,    'Dec': 31}

        month_index  = {'Jan' :1, 'Feb' :2,    'Mar' :3 ,'Apr' :4, 'May' :5, 'Jun' :6, \
                        'Jul' :7, 'Aug':8,   'Sep':9, 'Oct':10,  'Nov':11,  'Dec':12}     
            
        DJF = [10, 11, 0, 1]
        MAM = [1, 2, 3, 4]
        JJA = [4, 5, 6, 7]
        SON = [7, 8, 9, 10]   
            
        input_file   = dict(\
                             YEAR        = 0,\
                             ZEROTH_YEAR = 1979,\
                             START_YEAR  = begin_year  - 1979, \
                             END_YEAR    = finish_year - 1979,\
                             WINDOW_SIZE = 20,\
                             EVOLVE_WINDOW = 11,\
                             STARTING_DAYS_FOR_PREVIOUS_MONTH = 20,\
                             THRESHOLD   = 0.6,\
                             MONTH_INDEX = DJF,\
                             HEMISPHERE  = 'NH',\
                             FIELD_NAME  = 'A_N_daily',\
                             CRITICAL    =  65,\
                             source_NH   = '/data/pragallva/2023_repeat_ERA5/post_processing/parameters/',\
                             DESTINATION = '/data/pragallva/2023_repeat_ERA5/post_processing/blocking_info_extra_variable/', \
                             R           = 6378e3) 
        input_file['DESTINATION'] = input_file['DESTINATION']+'/DJF_Ac=%d/'%(input_file['CRITICAL'])
        h5saveload.make_sure_path_exists(input_file['DESTINATION'])
        
        for key in input_file.keys():
            globals()[key] = input_file[key]

        ############################# Load all the variables ################################## 
#         files_to_load = { 'A_N.hdf5'      : 'A_N_daily',   }        
#                           'U_N.hdf5'      : 'U_N_daily',   \
#                           'F1_N.hdf5'     : 'F1_N_daily',  \
#                           'F2_N.hdf5'     : 'F2_N_daily',  \
#                           'F3_N.hdf5'     : 'F3_N_daily',  \
#                           'Uref_N.hdf5'   : 'URB_N_daily', }

        logging_object0 = logruns.default_log(logfilename='loading', log_directory = './logs/') 
                
        A_N       = h5saveload.load_dict_from_hdf5(source_NH+'A_N.hdf5')
        A_N_daily = A_N['field']
        lat, lon  = A_N['lat'], A_N['lon']
        logging_object0.write('Extracted %s'%('A_N'))
        
        U_N_daily   = h5saveload.load_dict_from_hdf5(source_NH+'U_N.hdf5')   ['field'] 
        F2_N_daily  = h5saveload.load_dict_from_hdf5(source_NH+'F2a_N.hdf5')  ['field'] 
        F1_N_daily  = h5saveload.load_dict_from_hdf5(source_NH+'F1_N.hdf5')  ['field']
        F3_N_daily  = h5saveload.load_dict_from_hdf5(source_NH+'F3_N.hdf5')  ['field'] 
        URB_N_daily = h5saveload.load_dict_from_hdf5(source_NH+'Uref_N.hdf5')['field'] 
        
        F_N_daily      = (F1_N_daily + F2_N_daily + F3_N_daily)
        Uratio_N_daily = (U_N_daily/URB_N_daily)
        ############################# Load all the variables ################################## 

        ####################### lon extended and lat_mean #####################               
        lon_extended = np.append(lon, lon[:int(len(lon)/2)]+360)
        lat_mean     = lat.reshape(len(lat)//2,2).mean(axis=-1)
        #######################################################################
        
        for THIS_YEAR in range(START_YEAR,END_YEAR):
        
                logging_object = logruns.default_log(logfilename  = 'TRACK_test_%d'%(THIS_YEAR + ZEROTH_YEAR),  log_directory = './logs/')
                logging_object.write('----------------YEAR = %d-----------------------' %(THIS_YEAR + ZEROTH_YEAR))
                
                input_file['YEAR'] = THIS_YEAR 
                    
                FIELD = A_N_daily
                        
                for key in input_file.keys():
                    vars()[key] = input_file[key]
                
                if not os.path.exists(DESTINATION+'%d.hkl'%(THIS_YEAR+ZEROTH_YEAR)):
            
                        all_dates = return_window_dates(month_index=MONTH_INDEX, orig_year=YEAR, window_size=WINDOW_SIZE,
                                                        starting_day_for_previous_month=STARTING_DAYS_FOR_PREVIOUS_MONTH)[1]

                        logging_object.write('Total no of days being evaluated = %d' % len(all_dates))

                        save_data = {}

                        start = ti.time()
                        for START_DAY in range(len(all_dates)):
                            logging_object.write('Thats the start date %s, right?' % (all_dates[START_DAY][0]))
                            save_data, date = save_blocking_info(start_day=START_DAY, field=FIELD,
                                                                 critical=CRITICAL,fill=False, plot_boolean=False, window_size_days=WINDOW_SIZE,
                                                                 threshold=THRESHOLD,  year=YEAR, starting_day_for_previous_month=STARTING_DAYS_FOR_PREVIOUS_MONTH,
                                                                 month_index=MONTH_INDEX, lon=lon, lat=lat, evolve_window=EVOLVE_WINDOW, input_dict=save_data, \
                                                                 save_separate=True, DEST=DESTINATION)

                        end = ti.time()
                        logging_object.write('Time take for one year = %1.2f' % (end-start))

                        logging_object.write('----------Going to remove redundant blocks---------------------')
                        ### Define some functions to remove redundants #######


                        BLOCKS_IN_A_YEAR=remove_redundant_blocks(save_data) #### Somehow this does not manage to remove all redundants ####
                        logging_object.write('---------- Removed all redundant blocks -------------------')

                        #####  Save data ######

                        h5saveload.make_sure_path_exists(DESTINATION)
                        h5saveload.make_sure_path_exists(DESTINATION+'/input_file/')
                        hkl.dump(BLOCKS_IN_A_YEAR, DESTINATION+'%d.hkl' % (YEAR+ZEROTH_YEAR))
                        hkl.dump(input_file ,      DESTINATION+'/input_file/'+'%d.hkl' % (YEAR+ZEROTH_YEAR))
                        #######################

                        logging_object.write('Data written to disk')
                        logging_object.write('Great! Finished successfully')
                        
                else:
                                            
                        logging_object.write( DESTINATION+'%d.hkl exists already'%(YEAR+ZEROTH_YEAR) )
                                                
        #################################################################################################### 
        ############## Reload all data and remove any remaining duplicates #################################
        ####################################################################################################
        
        
        logging_object0.write('Reloading all data again') 
        for year in range(begin_year,finish_year):
                vars()["blocks%i"%year] = return_sorted_blocks(year, blocking_data = DESTINATION)
                
                
        logging_object0.write('Finding redundant data')         
        DELETE_INDEX={}
        for year in range(begin_year,finish_year):
            for keys_pair in itertools.combinations(eval("blocks%i"%year).keys(), 2):
                 DELETE_INDEX=compare_and_delete(keys_pair, eval("blocks%i"%year), DELETE_INDEX, year)
                    
        logging_object0.write('Deleting redundant data')            
        for key, values in (DELETE_INDEX.items()):
            year = values[0]   
            deletion_key=key
            list_index=values[1]

            all_keys=list(eval("blocks%i"%year)[deletion_key].keys())
            all_keys.remove('start_block_date')
            all_keys.remove('start_block_day')    

            for every_key in all_keys:
                eval("blocks%i"%year)[deletion_key][every_key].pop(list_index)  # <<<======= DELETION LINES FROM DICTIONARY  =======>>> 

            any_key='peak_block_date'
            if (len(eval("blocks%i"%year)[deletion_key]))==0:
                eval("blocks%i"%year).pop(key, None)         # <<<======= DELETION LINES FROM DICTIONARY  =======>>> 
                
            BLOCKS_IN_A_YEAR = eval("blocks%i"%year)
            hkl.dump(BLOCKS_IN_A_YEAR, DESTINATION+'%d.hkl' % (year))
        
        logging_object0.write('Great!! Saved all data. Today is %s. You had a good today!'%(cur_date))
        
        #################################################################################################### 
        ####################################################################################################
        ####################################################################################################        
        
        
        
                        
