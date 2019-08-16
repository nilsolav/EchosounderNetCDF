# this is the main test script for creating the NC interpretaion masks

# Region example
import datetime as dt
import numpy as np

# Super simple polygon defined per time. There are more fancy
# stuff that we could use, so please come up with suggestions 
time0 = dt.datetime(2019, 8, 16, 12, 10, 0, 1)
time = [time0 + dt.timedelta(seconds=x) for x in range(0, 10)]
# the depth contains the upper and lower bound
depth = [[10, 10, 10, 11, 12], [25, 26, 29, 40, 40]]

# NetCDF structure 

# \labels\labelID - vector with unique ID for the school
# \labels\starttime - vector with min(time) for each bounding box
# \labels\stoptime - vector with max(time) for each bounding box
# \labels\category - the acoustic category for this label

# \labels\label[1]\time - The time points for each label
# \labels\label[1]\depth - The depth range (start and stop) for each label
