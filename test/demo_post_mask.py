# -*- coding: utf-8 -*-
import h5py
import matplotlib.pyplot as plt
# from pathlib import Path
# import cftime

# Stuff to make the plots have nicely formated time axes
import datetime
import matplotlib.dates as mdates
import matplotlib.units as munits
import pdb
import os
import requests
import numpy as np
import glob

# baseUrl : local host of LSSS
baseUrl = 'http://localhost:8000'

# Name of the netcdf file
#direc = '//home//user//repos//echo-stuffs//'
#direc = 'D:\\DATA\\'
direc ='/mnt/d/DATA/'
#filenames = [direc +
#             'LSSS-label-versioning//' +
#             'S2016837//ACOUSTIC//LSSS//' +
#             'WORK//2016837-D20160427-T221032.nc']
dirname = [direc +
'S2018823_PEROS_3317/ACOUSTIC/LSSS/WORK/*.nc']
filenames=glob.glob(dirname[0])[100:-1]

def post_nc_lsss(filename, is_save_png=True, is_show=False):

    with h5py.File(filename, 'r') as f:
        # Open the group where the data is located
        interp = f['Interpretation/v1']
        
        # Get some variables and attributes
        t = interp['mask_times']
        d = interp['mask_depths']

        #d_units = 'm'  # str(interp['mask_depths'].attrs['units'], 'utf-8')
        #t_units = 's'  # str(interp['mask_times'].attrs['units'], 'utf-8')
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
        #time_fixer = 10000  # divide all times by this before using cftime.num2pydate()

        #cat_names = interp['region_category_names']
        #cat_prop = interp['region_category_proportions']
        #cat_ids = interp['region_category_ids']

        #for i, r in enumerate(cat_names):
        #    print('Region ' + str(cat_ids[i]) + ' has category '
        #          + '"' + cat_names[i] + '"'
        #          + ' with proportion ' + str(cat_prop[i]))
        post_masks(d, t)

def getFiletime(dt):
    # Convert from windows NTtime to python time
    microseconds = dt / 10
    seconds, microseconds = divmod(microseconds, 1000000)
    days, seconds = divmod(seconds, 86400)
    dtime = datetime.datetime(1601, 1, 1) + datetime.timedelta(
        days, seconds, microseconds)
    return dtime

# get masks
def post_masks(d, t):
    # Set zoom based on the time and depth from the data
    url = baseUrl + '/lsss/module/PelagicEchogramModule/zoom'
    min_time = min([min(r) for r in t])
    max_time = max([max(r) for r in t])
    min_dep = min([min(r) for r in d])
    max_dep = max([max(r) for r in d])
    zoom_req = [{"time":getFiletime(min_time).isoformat()+'Z',
                 "z":str(min_dep)}, {"time":getFiletime(max_time).isoformat()+'Z', "z":str(max_dep)}]
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

    # Convert raw data ping times to timestamps (Requires Python 3.7+ for datetime.isoformat())
    Tping_all = np.asarray([datetime.datetime.fromisoformat(tt[:-1]).timestamp() for tt in T])
    Pping_all = np.asarray(P)
    
    # Loop over schools/annotations
    for i, Rmask_all in enumerate(d):
        # Convert mask times from NT time to timestamps
        Tmask_all = [getFiletime(tt).timestamp() for tt in t[i]]
        
        # Restructure ranges to start-stop pairs by time
        Rmask_all_start = Rmask_all[0::2]
        Rmask_all_stop = Rmask_all[1::2]
        
        # since the format allow duplciate time for each depth range,
        # get the unique values for this school
        Tmask = np.unique(Tmask_all)
        
        # Find the time interval for this school
        Tping = []
        Tping.append(min(Tmask))
        Tping = Tping + (Tping_all[(Tping_all >= min(Tmask)) & (Tping_all <= max(Tmask))]).tolist()

        # Fill the json string
        json_str = []
        # Generate the mask string, loop over time
        for i, Tmask_i in enumerate(Tmask):
            min_max_vals = []
            # Loop over min max vals
            ind = Tmask_all==Tmask_i
            Rmask_i_start = Rmask_all_start[Tmask_all == Tmask_i]
            Rmask_i_stop = Rmask_all_stop[Tmask_all == Tmask_i]

            # Test if Rmask_i is empty, and if yes, then interpolate and
            # generate an interpolated Rmask_ii for this time step. this
            # needs to look at the pro and pre ceding ping to map out any
            # holes or gaps in the mask
            if len(Rmask_i_start)==0:
                print('Missing depth range for time step in the data. Must interpolate')
                pdb.set_trace()

            # Loop over the Rmask_i and add to json_str
            for j, Rmask_ii_start in enumerate(Rmask_i_start):
                min_max_vals.append({"min": Rmask_ii_start, "max": Rmask_i_stop[j]})

            # Add time and range to json string
            json_str.append({"time": datetime.datetime.fromtimestamp(Tmask_i).isoformat()+'Z',
                             "depthRanges": min_max_vals})
        
        # Post the school into LSSS
        if json_str:
            # print(*json_str, sep="\n")
            post('/lsss/module/PelagicEchogramModule/school-mask',
                 json = json_str)
        # Get the channel id
        # url2 = baseUrl + '/lsss/data/frequencies'
        # sounder_info = requests.get(url2).json()
        # freq = 38000.0
        # channel_id = [j for j,i in enumerate(sounder_info) if i==freq]
        # json_str2 = {'channels': [{'channelId': channel_id,
        #                           'categories': [{'id': 27, 'initials': 'Sand eel', 'assignment': 1}]}]}
        # region selection, få index
        # /lsss/regions/selection
        # loop
        # sett kategori
        
        # pdb.set_trace()
        # post('/lsss/module/InterpretationModule/scrutiny', json = json_str2)
        # Example http://localhost:8000/lsss/module/InterpretationModule/scrutiny data/pings?pingNumber=10&pingCount=100


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

# Run the shit
# Set the channel ID
# requests.post(baseUrl + '/lsss/data/frequency',json={'value':38000.0})
for filename in filenames:
    post_nc_lsss(filename)
