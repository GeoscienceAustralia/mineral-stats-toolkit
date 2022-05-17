# -*- coding: utf-8 -*-
"""
   Copyright 2022 Commonwealth of Australia (Geoscience Australia)

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

from toolkit.merge_models import MergeModels
from toolkit.seed_random_points import RandomPoints
import os
import numpy as np
import argparse


### start by defining some variables
# working directory with resistivity models, contains a folder for each 
# resistivity model
wd_inputs = r'C:\tmp\inputs'

if not os.path.exists(wd_inputs):
    os.mkdir(wd_inputs)

# define depths to interpolate model onto
depth_list = np.hstack([np.arange(0,60,2), 
                        np.arange(60,150,5)])*1e3

# number of repeats to run for random points
n_random_points = 20000 # number of random points to seed
                        # needs to be at least n_repeats * number of points 
                        # where n_repeats is the number of sets of random
                        # locations to seed and number of points is the number
                        # of deposits in the point dataset

# dictionary with details of ModEM inversions. 
# keys: basename of folder containing files (folder should be located within wd_resmodels)
# values: list containing:
#    - model file name (*.rho file)
#    - centre point of model [x,y] in projected coordinates
#    - proj4 string for projection
modem_details_dict = {
                'Robertson2020':['Delamerian_AusLAMP_cropped.rho',[366832.99996643455, 6428315.2838392835],
                                '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                'Kirkby2020':['AusLAMP_NSW_Vic.rho',[3.789018683870512177e+05,6.158864711905587465e+06],
                                     '+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
               }

# first, merge models
resmodel = MergeModels(wd_inputs,
                       savepath = wd_inputs,
                       cellsize=0.1,
                       station_buffer=0.7,
                       stationxy_fn_basename = 'Station_XY.txt',
                       modem_details_dict=modem_details_dict,
                       savefilename_prefix='resmodels',
                       depth_list=depth_list)
resmodel.create_netcdf()


# second, create some random points covering same area as resistivity models
rp = RandomPoints(wd_inputs, 
                  modem_details_dict=modem_details_dict,
                  savepath=wd_inputs,
                  n_points=n_random_points,
                  chunksize=10000000)
rp.seed_random_points()