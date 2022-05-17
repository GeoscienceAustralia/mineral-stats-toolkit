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

Script to run an analysis from start to finish.
See script create_inputs_from* to generate inputs from your format

"""

from toolkit.distance_to_contour import DistanceToContour
from toolkit.cdf import CDF
import os
import argparse

parser = argparse.ArgumentParser(description='Command line input.')
parser.add_argument("--mst", action="store", help="Path to the Mineral-Stats-Toolkit installation directory.")
rv = parser.parse_args()

### start by defining some variables
# working directory 
wd_inputs = r'C:\tmp\inputs'
wd_outputs = r'C:\tmp\outputs'


if not os.path.exists(wd_outputs):
    os.mkdir(wd_outputs)

# define paths to input files,
# points file (containing deposit location and size etc)
point_file = os.path.join(wd_inputs,'Orogenic_gold.txt')
# netcdf file with geoscience dataset to analyse
netcdf_fn = os.path.join(wd_inputs,'resmodels.nc')
# station locations
station_locations_fn = os.path.join(wd_inputs,'resmodels_sloc.txt')
# random locations
random_locations_fn = os.path.join(wd_inputs,'RandomPoints_20000points.npy')

# define target value (in this case 100 ohm-m)
target_value = 100

# number of repeats to run for random points
n_repeats = 100

# first, calculate distances from contour        
distanceObj = DistanceToContour(netcdf_fn,
                                station_locations_fn,
                                point_file, 
                                target_value,
                                random_locations_fn = random_locations_fn,
                                savepath = wd_outputs,
                                point_propertynames=['Au','Cu','Age_Ga'],
                                n_repeats=n_repeats)
distanceObj.setup_and_run_analysis()
