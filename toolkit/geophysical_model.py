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

import os
import numpy as np
from toolkit.modem import ModEM
from scipy.spatial import ConvexHull
from shapely.geometry import MultiPoint,Polygon,MultiPolygon

class GeophysicalModel():
    def __init__(self, wd, **kwargs):
        """
        
        class to contain details on geophysical models and station locations
        used in geophysical analysis

        Parameters
        ----------
        wd : TYPE
            DESCRIPTION.
            
            
        Optional parameters
        -------------------
        station_buffer : buffer around stations to save resistivity values
        model_type : type of models to combine into netcdf. Default is ModEM.
            Currently, only ModEM is supported. However, there is the
            possibility to include other data types. Need to build a class to
            sit around other data type, equivalent to the ModEM class used 
            here. The class needs to have attributes gcx, gcy, gcz (grid 
            centres), slon, slat (station lat/lon) and values (physical 
            property, 3D array with values corresponding to the gcx, gcy and 
            gcz locations). Also need to add a function within the 
            GeophysicalModel class to read in the model
        modem_details_dict : dict
            Dictionary containing parameters of ModEM inversions. Keys are the
            subfolders within wd in which the ModEM inversion files are located.
            Values are lists containing three values:
                modem_rho_file : resistivity model in .rho format as output by
                                 ModEM
                centre_point : list containing the x, y coordinate of the
                               zero point in the model grid, in real world
                               coordinates defined by proj_str
                proj_str : Proj4 string (as per spatialreference.info) to
                           define the projection that model is in.
            example :
                modem_details_dict = {'model1':['modem_rho_file1', 
                                                [650000,4320000], 
                                                '+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                                      'model2':['modem_rho_file2', 
                                                [0,5220000], 
                                                '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                                      }

        stationxy_fn_basename : if ModEM format is used, need to provide
            a column based file with their locations in local grid coordinates.
            this argument provides the name of this file.
        propertyname : name of the physical property in the geophysical model
            (e.g. resistivity)


        Returns
        -------
        None.

        """
        
        self.wd = wd
        self.station_buffer = 0.7
        
        self.model_type = 'ModEM'
        self.modem_details_dict = None
        self.stationxy_fn_basename = 'Station_XY.txt'
        self.propertyname = None
        
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self,key,kwargs[key])

    def read_models(self):
        if self.model_type == 'ModEM':
            self.propertyname = 'resistivity'
            self.read_modem()

    def read_modem(self):
        """
        
        Read ModEM models and into a dictionary stored as an attribute

        Returns
        -------
        None.

        """
        self.model_dict = {}

        from toolkit.modem import ModEM
        
        
        for key in self.modem_details_dict.keys():
            # print("reading",key)
                
            
            model_fn, centre, proj_str = self.modem_details_dict[key]
            
            modObj = ModEM(os.path.join(self.wd,key,model_fn),
                           proj_str=proj_str,
                           stationxy_fn = os.path.join(self.wd, key,
                                                       self.stationxy_fn_basename),
                           epsg=0,
                           rpad=1,
                           centre=centre)
            
            
            modObj.mask_station_distance(self.station_buffer)
                
            self.model_dict[key] = modObj
            
            

        
    def compile_station_locations(self):
        """
        Compile station x, y (longitude, latitude) into a list.
        Also store 
        - a station multipoint (for all stations)
        - convex hull around the stations for each model (multipolygon)
        - min/max lat/lon station locations for each model (dict with keys as
                                                            for model_dict)

        Returns
        -------
        None.

        """
        self.station_ll = None
        station_polylist = []
        self.bounds_dict = {}
        
        for key in self.model_dict.keys():
            modObj = self.model_dict[key]
            # print("compiling stations for",key)
            
            # list of station latitude and longitude
            sll_list = np.stack([modObj.slon,modObj.slat],axis=1)
            
            
            # append to existing station lat/lon
            if self.station_ll is None:
                self.station_ll = sll_list
            else:
                self.station_ll = np.concatenate([self.station_ll,sll_list])
                
            # calculate and store convex hull around stations
            hull =ConvexHull(sll_list)
            station_polylist.append(Polygon(sll_list[hull.vertices]))
            
            # calculate and store min/max lat/lon of stations
            minLon = modObj.slon.min()
            minLat = modObj.slat.min()
            maxLon = modObj.slon.max()
            maxLat = modObj.slat.max()
            self.bounds_dict[key] = [minLon, minLat, maxLon, maxLat]
            
        self.station_multipoint = MultiPoint(self.station_ll)
        self.station_multipolygon = MultiPolygon(station_polylist)
        self.station_polylist = station_polylist