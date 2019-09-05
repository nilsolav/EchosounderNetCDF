# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:44:06 2019

@author: gavinj
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt

#Name of the netcdf file
filename = 'demo_mask_vlen.nc'

f = h5py.File(filename, 'r')

#Open the group where the data is located
interp = f['Interpretation/v1']

# Get some variables and attributes
t = np.array(interp['mask_times'])
d = np.array(interp['mask_depths'])

d_units = str(interp['mask_depths'].attrs['units'], 'utf-8')

t_units = str(interp['mask_times'].attrs['units'], 'utf-8')
t_calendar = str(interp['mask_times'].attrs['calendar'], 'utf-8')

c = interp['sound_speed'][()]
c_units = str(interp['sound_speed'].attrs['units'], 'utf-8')

bb_upper = np.array(interp['min_depth'])
bb_lower = np.array(interp['max_depth'])
bb_left = np.array(interp['start_time'])
bb_right = np.array(interp['end_time'])

region_id = np.array(interp['id'])
region_name = np.array(interp['name'])
r_type = np.array(interp['region_type'])
r_type_enum = h5py.check_dtype(enum=interp['region_type'].dtype)
# Convert region types into a text version
r_type_enum = dict(map(reversed, r_type_enum.items()))
r_type_name = [r_type_enum[i] for i in r_type]

#close dataset
f.close()

#Plot the power of beam
plt.figure(1)
plt.clf()

# plot masks
for i, r in enumerate(d):
    for time, ranges in zip(t[i], r):
        print(time)
        for start, stop in zip(ranges[0::2], ranges[1::2]):
            plt.plot(np.array([time, time]), np.array([start, stop]), linewidth=4, color='k')

# plot bounding boxes
for i, junk in enumerate(bb_upper):
    plt.plot(np.array([bb_left[i], bb_right[i], bb_right[i], bb_left[i], bb_left[i]]), 
             np.array([bb_lower[i], bb_lower[i], bb_upper[i], bb_upper[i], bb_lower[i]]),
             color=(0.5, 0.5, 0.5))
    plt.text(bb_left[i], bb_upper[i], 'ID: ' + str(region_id[i]) 
            + ' (' + r_type_name[i] + ')', 
             bbox=dict(facecolor=(0.5, 0.5, 0.5), alpha=0.5))
plt.title('Using c= ' + str(c) + ' ' + c_units)
ax = plt.gca()
ax.invert_yaxis()

plt.xlabel('Time\n(' + t_units + ')')
plt.ylabel('Depth (' + d_units + ')')

plt.savefig('demo_mask_vlen.png', bbox_inches='tight')
