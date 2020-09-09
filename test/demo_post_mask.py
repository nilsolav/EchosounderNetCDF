# -*- coding: utf-8 -*-
"""
1) This work is based on the initial work that was created on
   Tue Aug 20 11:44:06 2019 by gavinj that generates masks
of 3 rectangles
2) Debugging were furthered by nils olav in july 2020
3) This work is then integrated with functions to post via
   LSSS API by yi liu in july 2020
"""

# Gavin/Nillav - Read netcdf file
import h5py
import matplotlib.pyplot as plt
# from pathlib import Path
# import cftime

# Stuff to make the plots have nicely formated time axes
import datetime
import matplotlib.dates as mdates
import matplotlib.units as munits
import pdb
# Yi - posting and integration
import os
import requests
import numpy as np
# import scipy # converts polygon to ping-based formats
# import math
#from datetime import datetime

converter = mdates.ConciseDateConverter()
munits.registry[np.datetime64] = converter
munits.registry[datetime.date] = converter
munits.registry[datetime.datetime] = converter

# Global variables, can be changed in docker

# baseUrl : local host of LSSS
baseUrl = 'http://localhost:8000'

# Name of the netcdf file - experimental mask data in NetCDF4 format
#direc = '//home//user//repos//echo-stuffs//'
#direc = 'D:\\DATA\\'
direc ='/mnt/d/DATA/'

filenames = [direc +
             'LSSS-label-versioning//' +
             'S2016837//ACOUSTIC//LSSS//' +
             'WORK//2016837-D20160427-T221032.nc']

''' 
Function areas
nc_reader : function that reads masks along with the attributes from nc file
post      : function that post the data (mask) onto LSSS given path and parameters
'''


def nc_reader(filename, is_save_png=True, is_show=False):

    with h5py.File(filename, 'r') as f:

        # Open the group where the data is located
        interp = f['Interpretation/v1']

        # Get some variables and attributes
        t = interp['mask_times']
        d = interp['mask_depths']
        d_units = 'm'  # str(interp['mask_depths'].attrs['units'], 'utf-8')
        t_units = 's'  # str(interp['mask_times'].attrs['units'], 'utf-8')
        #t_calendar = str(interp['mask_times'].attrs['calendar'], 'utf-8')

        #c = interp['sound_speed'][()]
        #c_units = str(interp['sound_speed'].attrs['units'], 'utf-8')

        #bb_upper = interp['min_depth']
        #bb_lower = interp['max_depth']
        #bb_left = interp['start_time']
        #bb_right = interp['end_time']
        #region_id = interp['region_id']
        # region_name = interp['region_name']
        #r_type = interp['region_type']
        #r_type_enum = h5py.check_dtype(enum=interp['region_type'].dtype)

        # Convert region types into a text version
        #r_type_enum = dict(map(reversed, r_type_enum.items()))
        #r_type_name = [r_type_enum[i] for i in r_type]

        # convert time variables into the form that matplotlib wants
        # current example .nc files actually have timestamps in 100 nanseconds since 1601.
        # But give the time units as milliseconds since 1601. Sort that out...
        time_fixer = 10000  # divide all times by this before using cftime.num2pydate()

        #cat_names = interp['region_category_names']
        #cat_prop = interp['region_category_proportions']
        #cat_ids = interp['region_category_ids']

        #for i, r in enumerate(cat_names):
        #    print('Region ' + str(cat_ids[i]) + ' has category '
        #          + '"' + cat_names[i] + '"'
        #          + ' with proportion ' + str(cat_prop[i]))

        if is_show or is_save_png:
            #plt.figure()
            #plt.clf()
            get_masks(d, t, d_units, t_units, time_fixer, handle=plt)
            #get_region(region_id, bb_upper, bb_lower, bb_left, bb_right,
            #           r_type_name, time_fixer, handle=plt)

            # Plot the power of beam
            #plt.title('Using c= ' + str(c) + ' ' + c_units)
            #ax = plt.gca()
            #ax.invert_yaxis()

            #plt.xlabel('Time\n(' + t_units + ')')
            #plt.ylabel('Depth (' + d_units + ')')

            #if is_save_png:

            #    plt.savefig(os.path.splitext(filename)[0] +
            #                '.png', bbox_inches='tight')

            #if is_show:
            #    plt.show()
        else:
            get_masks(d, t, d_units, t_units, time_fixer, handle=None)
            #get_region(region_id, bb_upper, bb_lower, bb_left,
            #           bb_right, r_type_name, time_fixer, handle=None)


def getFiletime(dt):
    # Convert from windows NTtime to python time
    microseconds = dt / 10
    seconds, microseconds = divmod(microseconds, 1000000)
    days, seconds = divmod(seconds, 86400)
    dtime = datetime.datetime(1601, 1, 1) + datetime.timedelta(
        days, seconds, microseconds)
    return dtime

# get masks
def get_masks(d, t, d_units, t_units, time_fixer, handle=None):
    # Set zoom based on the time and depth from the data
    url = baseUrl + '/lsss/module/PelagicEchogramModule/zoom'
    min_time = min([min(r) for r in t])
    max_time = max([max(r) for r in t])
    min_dep = min([min(r) for r in d])
    max_dep = max([max(r) for r in d])

    zoom_req = [{"time":getFiletime(min_time).isoformat()+'Z', "z":str(min_dep)}, {"time":getFiletime(max_time).isoformat()+'Z', "z":str(max_dep)}]
    post('/lsss/module/PelagicEchogramModule/zoom', json=zoom_req)

    # Get time and ping array from visible view on echogram
    url = baseUrl + '/lsss/module/PelagicEchogramModule/zoom'
    response = requests.get(url).json()
    t0 = response[0]['pingNumber']
    t1 = response[1]['pingNumber']
    tint = t1-t0+1
    # Get the ping/time array for interpolation from region to LSSS shapes
    # Example http://localhost:8000/lsss/data/pings?pingNumber=10&pingCount=100
    url = baseUrl + '/lsss/data/pings?pingNumber=' + str(t0)+'&pingCount=' + str(tint)
    # Get the list
    ping_time = requests.get(url).json()
    # Get time and ping lists (for interpolating empty pings)
    T = [nilz['time'] for nilz in ping_time]
    P = [nilz['pingNumber'] for nilz in ping_time]
    
    # Loop over schools/annotations
    for i, Rmask_all in enumerate(d):
        # Convert mask times from NT time to timestamps
        Tmask_all = [getFiletime(tt).timestamp() for tt in t[i]]
        #pdb.set_trace()
        # Restructure ranges to start-stop pairs by time
        d_start = Rmask_all[0::2]
        d_stop = Rmask_all[1::2]
        
        # Convert r values array to list
        pointsmin = min(Rmask_all)
        pointsmax = max(Rmask_all)

        # since the format may have multiple times, get the unique values for this school
        Tmask = np.unique(Tmask_all)

        # Get the times where there are multiple range intervals for each time step
        
        diffs = ((np.diff(Tmask_all) == 0)*int(1))
        indices = np.append(1,np.cumsum(diffs))

        diffinds = np.where(diffs)[0]
        #nilz = r[diffinds]

        # Convert raw data ping times to timestamps (Requires Python 3.7+ for datetime.isoformat())
        Tping_all = np.asarray([datetime.datetime.fromisoformat(tt[:-1]).timestamp() for tt in T])
        
        Tping = []
        # Add time start
        Tping.append(min(Tmask))
        # Add bisection
        Tping = Tping + (Tping_all[(Tping_all >= min(Tmask)) & (Tping_all <= max(Tmask))]).tolist()
        # Add time end
        # Tping.append(max(Tmask))

        # Interpolate min and max
        #vals_min = np.interp(Tping, Tmask, pointsmin)
        #vals_max = np.interp(Tping, Tmask, pointsmax)

        json_str = []
        # Yi: This part needs to be changed to add in missing time steps.

        #for time, ranges in zip(t[i], r):
        #    min_max_vals = []
        #    # Ranges
        #    # pdb.set_trace()
        #    for start, stop in zip(ranges[0::2], ranges[1::2]):
        #        min_max_vals.append({"min": float(start), "max": float(stop)})

        #        if handle:
        #            handle.plot([time, time], [start, stop],
        #                        linewidth=4, color='k')

        # Append original + bisected (interpolated) time
        #for xx, vv in enumerate(bisected):
        #    min_max_vals = []
        #    min_max_vals.append({"min": vals_min[xx], "max": vals_max[xx]})
        #    json_str.append({"time": datetime.fromtimestamp(vv).isoformat()+'Z',
        #                     "depthRanges": min_max_vals})

        # Generate the mask string
        for xx, vv in enumerate(Tmask_all):
            min_max_vals = []
            min_max_vals.append({"min": d_start[xx], "max": d_stop[xx]})
            json_str.append({"time": datetime.datetime.fromtimestamp(vv).isoformat()+'Z',
                             "depthRanges": min_max_vals})
        # Post the school into LSSS
        if json_str:
            #print(*json_str, sep="\n")
            post('/lsss/module/PelagicEchogramModule/school-mask',
                        json = json_str)

# get region
def get_region(region_id, bb_upper, bb_lower, bb_left, bb_right, r_type_name, time_fixer, handle = None):

    for i, junk in enumerate(bb_upper):

        if handle:
            plt.plot([bb_left[i], bb_right[i], bb_right[i], bb_left[i], bb_left[i]],
                     [bb_lower[i], bb_lower[i], bb_upper[i], bb_upper[i], bb_lower[i]],
                     color=(0.5, 0.5, 0.5))

            plt.text(bb_left[i], bb_upper[i], 'ID: ' + str(region_id[i])
                    + ' (' + r_type_name[i] + ')',
                     bbox=dict(facecolor=(0.5, 0.5, 0.5), alpha=0.5))


def post(path, params=None, json=None, data=None):
    '''
    This is the basic function to post a mask/masks to LSSS via the API. Following
    is the parameters
    - path   : directory/category of a mask/masks to post to
    - params : parameters specifying the post process, default is None
    - json   : the mask is in jason format, in this case, is a list of dictionary for example
             [ {'pingNumber': 54001, 'depthRanges': [{'min': 30.6, 'max': 34.1}]},
               {'pingNumber': 54002, 'depthRanges': [{'min': 30.6, 'max': 37.5}]},
               {'pingNumber': 54003, 'depthRanges': [{'min': 30.6, 'max': 39.9}]},
                ...
             ]
             Note that
                * you can also post multiple "min"/"max" ranges in case the mask is not a rectangle
                * 'pingNumber' should be continuous integer, otherwise, the posted results will be multiple polygons
                * 'depthRange' can be float value or integer

    - data   : data (usually a ping) to post to the LSSS GUI, which is slow thus not recommended when the data is big.
               the same situation for the request command.

    '''

    # Make sure the LSSS server is on then the url points to the target directory of masks/areas/schools/layers
    url = baseUrl + path

    # Connect to the server and post the content
    response = requests.post(url, params=params, json=json, data=data)


    # Check the feedback from LSSS API
    if response.status_code == 200:
        return response.json()
    if response.status_code == 204:
        return None
    raise ValueError(url + ' returned status code ' + str(response.status_code) + ': ' + response.text)

if __name__ == "__main__":
    # Batch check all nc files
    # for i, f in enumerate(filenames):
    #     # print(i, f, '\n')
    #     nc_reader(f, is_save=False, is_show=False)

    # Check individual nc file
    print('Files: ', filenames)

    nc_reader(filenames[0], is_save_png=False, is_show=True)

 
