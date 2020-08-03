# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:44:06 2019

@author: gavinj
"""

import h5py
import matplotlib.pyplot as plt
from pathlib import Path
import cftime
import numpy as np
import os

# Stuff to make the plots have nicely formated time axes
import datetime
import matplotlib.dates as mdates
import matplotlib.units as munits
converter = mdates.ConciseDateConverter()
munits.registry[np.datetime64] = converter
munits.registry[datetime.date] = converter
munits.registry[datetime.datetime] = converter

#Name of the netcdf file
filename = 'demo_mask.nc'
baseDir = 'F:\\Documents\\local\\netcdf formats\\LSSS-label-versioning'
cruises = ['S2005114', 'S2006101', 'S2007205', 'S2009107', 'S2010205', 'S2012817',
           'S2012837', 'S2014119', 'S2016837', 'S2017838', 'S2017843']

# build up a list of nc files to process
ncfiles = []
for cruise in cruises:
    files = Path(baseDir, cruise, 'ACOUSTIC', 'LSSS', 'WORK').glob('*.nc')
    for file in files:
        ncfiles.append(file)

# for each .nc file, make an plot of the regions
for p in ncfiles:
    
    with h5py.File(p, 'r') as f:
    
        #Open the group where the data is located
        interp = f['Interpretation/v1']
        
        # Get some variables and attributes
        t = interp['mask_times']
        d = interp['mask_depths']
        
        d_units = str(interp['mask_depths'].attrs['units'], 'utf-8')
        
        t_units = str(interp['mask_times'].attrs['units'], 'utf-8')
        t_calendar = str(interp['mask_times'].attrs['calendar'], 'utf-8')
        
        c = interp['sound_speed'][()]
        c_units = str(interp['sound_speed'].attrs['units'], 'utf-8')
        
        bb_upper = interp['min_depth']
        bb_lower = interp['max_depth']
        bb_left = interp['start_time']
        bb_right = interp['end_time']
        
        region_id = interp['region_id']
        region_name = interp['region_name']
        r_type = interp['region_type']
        r_type_enum = h5py.check_dtype(enum=interp['region_type'].dtype)
        # Convert region types into a text version
        r_type_enum = dict(map(reversed, r_type_enum.items()))
        r_type_name = [r_type_enum[i] for i in r_type]
        
        # convert time variables into the form that matplotlib wants
        # current example .nc files actually have timestamps in 100 nanseconds since 1601.
        # But give the time units as milliseconds since 1601. Sort that out...
        time_fixer = 10000 # divide all times by this before using cftime.num2pydate()
        
        cat_names = interp['region_category_names']
        cat_prop = interp['region_category_proportions']
        cat_ids = interp['region_category_ids']
        
        for i, r in enumerate(cat_names):
            print('Region ' + str(cat_ids[i]) + ' has category '
                  + '"' + cat_names[i] + '"'
                  + ' with proportion ' + str(cat_prop[i]))
        
        #Plot the power of beam
        plt.figure(1)
        plt.clf()
        
        # plot masks
        for i, r in enumerate(d):
            for time, ranges in zip(t[i], r):
                time = cftime.num2pydate(time/time_fixer, t_units)
                for start, stop in zip(ranges[0::2], ranges[1::2]):
                    plt.plot([time, time], [start, stop], linewidth=4, color='k')
        
        # plot bounding boxes
        for i, junk in enumerate(bb_upper):
            left = cftime.num2pydate(bb_left[i]/time_fixer, t_units)
            right = cftime.num2pydate(bb_right[i]/time_fixer, t_units)
            plt.plot([left, right, right, left, left], 
                     [bb_lower[i], bb_lower[i], bb_upper[i], bb_upper[i], bb_lower[i]],
                     color=(0.5, 0.5, 0.5))
            plt.text(left, bb_upper[i], 'ID: ' + str(region_id[i]) 
                    + ' (' + r_type_name[i] + ')', 
                     bbox=dict(facecolor=(0.5, 0.5, 0.5), alpha=0.5))
            
        plt.title('Using c= ' + str(c) + ' ' + c_units)
        ax = plt.gca()
        ax.invert_yaxis()
        
        plt.xlabel('Time')
        plt.ylabel('Depth (' + d_units + ')')
        
        plt.savefig(p.with_suffix('.png'), bbox_inches='tight', dpi=300)
