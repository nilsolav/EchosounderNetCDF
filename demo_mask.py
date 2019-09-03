# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:44:06 2019

@author: gavinj
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt

#Name of the netcdf file
filename = 'demo_mask.nc'

#open the file
dataset = Dataset(filename)

#Open the group where the data is located
interp = dataset.groups['Interpretation']

# Get some variables and attributes
t = interp.variables['mask_times'][:]
d = interp.variables['mask_depths'][:]

c = interp.variables['sound_speed'][:]
c_units = interp.variables['sound_speed'].units

t_units = interp.variables['mask_times'].units
d_units = interp.variables['mask_depths'].units

bb_upper = interp.variables['min_depth'][:]
bb_lower = interp.variables['max_depth'][:]
bb_left = interp.variables['start_time'][:]
bb_right = interp.variables['end_time'][:]

region_id = interp.variables['id'][:]

#close dataset
dataset.close()

#Plot the power of beam
plt.figure(1)
plt.clf()

# plot masks
for i, r in enumerate(d):
    plt.plot(np.array([t[i], t[i]]), np.array([r[0::2], r[1::2]]), linewidth=4, color='k')

# plot bounding boxes
for i, junk in enumerate(bb_upper):
    plt.plot(np.array([bb_left[i], bb_right[i], bb_right[i], bb_left[i], bb_left[i]]), 
             np.array([bb_lower[i], bb_lower[i], bb_upper[i], bb_upper[i], bb_lower[i]]),
             color=(0.5, 0.5, 0.5))
    plt.text(bb_left[i], bb_upper[i], 'ID: ' + str(region_id[i]), 
             bbox=dict(facecolor=(0.5, 0.5, 0.5), alpha=0.5))
plt.title('Using c= ' + str(c) + ' ' + c_units)
ax = plt.gca()
ax.invert_yaxis()

plt.xlabel('Time\n(' + t_units + ')')
plt.ylabel('Depth (' + d_units + ')')

plt.savefig('demo_mask.png', bbox_inches='tight')
