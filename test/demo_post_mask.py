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
direc ='/home/user/repos/echo-stuffs/ncs'
#filenames = [direc +
#             'LSSS-label-versioning//' +
#             'S2016837//ACOUSTIC//LSSS//' +
#             'WORK//2016837-D20160427-T221032.nc']
dirname = [direc + '/*.nc']
filenames=sorted(glob.glob(dirname[0]))

def post_nc_lsss(filename, is_save_png=True, is_show=False):

    print(filename)
    with h5py.File(filename, 'r') as f:
        # Check if we have the times and depths
        if(('Interpretation/v1/mask_times' in f and 'Interpretation/v1/mask_depths' in f)):

            # Get some variables and attributes
            t = f['Interpretation/v1/mask_times']
            d = f['Interpretation/v1/mask_depths']

            # Post it
            post_masks(d, t)

def getFiletime(dt):
    # Convert from Windows NT FILETIME (nanoseconds from 1.1.1601) to Python date
    return (datetime.datetime(1601, 1, 1) + datetime.timedelta(microseconds=(dt / 10)))

# get masks
def post_masks(d, t):
    # Loop over all schools/annotations
    Tmask_all = np.array([])
    Rmask_all_start = np.array([])
    Rmask_all_stop = np.array([])

    for i, Rmask_all in enumerate(d):
        # Convert mask times from NT time to timestamps
        Tmask_all = np.concatenate([Tmask_all, [getFiletime(tt).timestamp() for tt in t[i]]])
        
        # Restructure ranges to start-stop pairs by time
        Rmask_all_start = np.concatenate([Rmask_all_start, Rmask_all[0::2]])
        Rmask_all_stop = np.concatenate([Rmask_all_stop, Rmask_all[1::2]])
        

    # since the format allow duplicate time for each depth range,
    # get the unique values for this school
    Tmask = np.unique(Tmask_all)
    
    # Fill the json string
    json_str = []
    # Generate the mask string, loop over time
    for i, Tmask_i in enumerate(Tmask):
        min_max_vals = []
        # Loop over min max vals
        ind = Tmask_all==Tmask_i
        Rmask_i_start = Rmask_all_start[ind]
        Rmask_i_stop = Rmask_all_stop[ind]

        # Test if Rmask_i is empty, and if yes, then interpolate and
        # generate an interpolated Rmask_ii for this time step. this
        # needs to look at the pro and pre ceding ping to map out any
        # holes or gaps in the mask
        if len(Rmask_i_start)==0:
            print('Missing depth range for time step in the data. Must interpolate')
            pdb.set_trace()

        # combine max min
        combined_masks = np.array([(x, y) for x,y in zip(Rmask_i_start, Rmask_i_stop)])

        # remove duplicates
        cleaned_masks = np.unique(np.sort(combined_masks, axis=1).view(','.join([combined_masks.dtype.char]*2))).view(combined_masks.dtype).reshape(-1, 2)

        # Loop over the Rmask_i and add to json_str
        for min, max in cleaned_masks:
            min_max_vals.append({"min": min, "max": max})

        # Add time and range to json string
        json_str.append({"time": datetime.datetime.fromtimestamp(Tmask_i).isoformat()+'Z',
                         "depthRanges": min_max_vals})
    
    # Post the school into LSSS
    post('/lsss/module/PelagicEchogramModule/school-mask',
         json = json_str, errorpass = True)
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


def post(path, params=None, json=None, data=None, errorpass=False):
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
        print("PASS:")
        print(*json, sep="\n")
        return response.json()
    if response.status_code == 204:
        return None

    print("FAILS:")
    print(*json, sep="\n")
    if errorpass == False:
        raise ValueError(url + ' returned status code ' + str(response.status_code) + ': ' + response.text)
    else:
        return None

if __name__ == "__main__":
    # Batch check all nc files
    # for i, f in enumerate(filenames):
    #     # print(i, f, '\n')
    #     nc_reader(f, is_save=False, is_show=False)

    # Check individual nc file
    print('Files: ', filenames)
    print('File: ', filenames[0])
    post_nc_lsss(filenames[0], is_save_png=False, is_show=True)

# Run the shit
# Set the channel ID
# requests.post(baseUrl + '/lsss/data/frequency',json={'value':38000.0})
for filename in filenames:
    post_nc_lsss(filename)
