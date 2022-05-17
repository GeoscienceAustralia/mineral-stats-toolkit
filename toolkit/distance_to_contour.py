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

from netCDF4 import Dataset
from toolkit.distance_to_contour_tools import epsg_project, filter_points, compute_contour,\
    compute_distance_to_target_contour
import numpy as np
from shapely.geometry import Point
import os
import time

class DistanceToContour():

    def __init__(self, netcdf_file, station_locations_file, 
                 point_file, target_value, **kwargs):
        """
        
        Class to manage calculation of the distance of point locations (e.g.
        mineral deposits) to a contour in a 3D gridded geoscientific dataset 
        (e.g. geophysical model). 
        
        Parameters
        ----------
        netcdf_file : str
            Full path to netcdf file containing geoscientific dataset to analyse.
            values in the file are longitude (1d array, shape l), latitude (1d
            array, shape m), depth (1d array, shape n), and values (l x m x n array)
        station_locations_file : str
            Full path to text file with station locations. File should be a
            space delimited column-based file with no header and columns lon, lat
        point_file : str
            Full path to point file containing points that will be analysed in
            terms of their distance to a contour in the geoscientifc dataset
            provided by netcdf_file. File needs to be space delimited with no
            header. First 2 columns are latitude (column 0) and longitude 
            (column 1). Can have any number of additional columns if desired
            and these will be filtered and saved into point_property_file. For
            example, point_properties could include metal tonnage or age data
            for mineral deposits
        target_value : float or int
            target value in geoscientific dataset to generate contours
            
        Optional Parameters
        ----------
        projection : str
            proj4 string for the projection to carry out the analysis in. The
            default is equal earth projection ('+proj=eqearth')
        savepath : str
            path to save results to
        random_locations_fn : filename for random locations (faster to pre-load
            than to calculate these on the fly)
        propertyname : name of property to use in the netcdf file
        n_repeats : number of realisations to calculate distances for in the 
            random locations
        outside_polygonfile : shapefile containing polygon/s within which to
            exclude points from analysis
        within_polygonfile : shapefile containing polygon/s outside which to
            exclude points from analysis
        station_buffer : resistivity values are only saved if they are located
            within station_buffer (in degrees) of an MT station in the 
            inversion. 
        prefix : prefix for the filenames to save the distances and associated
            data arrays to
        point_propertynames : names of additional properties in the point_data
            file. At this stage it is hardcoded to require 3 point_property
            columns
        Returns
        -------
        None.
        
        
        Example
        -------
        DistanceObj = DistanceToContour(netcdf_file, station_locations_file, 
                                        point_file, target_value)
        DistanceObj.setup_and_run_analysis()
        
        
        """
        self.netcdf_file = netcdf_file # netcdf file containing
        self.point_file = point_file
        self.target_value = target_value
        self.station_locations_file = station_locations_file

        
        self.projection = '+proj=eqearth'
        self.savepath = '.' 
        self.random_locations_fn = None 

        self.propertyname = 'resistivity'
        self.n_repeats = 100
        self.outside_polygonfile = None
        self.within_polygonfile = None
        self.station_buffer = 0.7 # station buffer in degrees latitude/longitude
        self.prefix = '%s_%1i'%(os.path.basename(point_file)[:-4],target_value)
        self.point_propertynames = None
        

        if self.propertyname == 'resistivity':
            self.prefix += 'ohmm'
        
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self,key,kwargs[key])
        
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)
            
        
        self.read_inputs()
        
        if self.point_propertynames is not None:
            self.point_propertynames = ['f%1i'%i for i in \
                            range(len(self.point_data[0])-2)]
                
            
    def setup_and_run_analysis(self):
        """
        
        Set up and run analysis
        Returns
        -------
        None.
        """
        
        self.read_inputs()
        self.prepare_data()
        self.run_analysis()


    def read_inputs(self):
        """
        
        Read latitude, longitude, and values from the netcdf file (defined by
        the netcdf_file attribute), point data (point_file attribute), 
        station locations (station_locations_file), and pre-calculated random
        locations (random_locations_file).
        
        """
        
        # read geophysical model/dataset
        ds = Dataset(self.netcdf_file)
        self.gclon, self.gclat = np.meshgrid(ds['longitude'][:],
                                             ds['latitude'][:])
        
        self.gcx, self.gcy = epsg_project(self.gclon,self.gclat,
                                          4326,0,
                                          proj_str=self.projection)
        self.gcz = ds['depth'][:]
        
        self.values = ds[self.propertyname][:]
        
        # read station locations
        self.station_ll = np.loadtxt(self.station_locations_file)
        
        # load points
        self.point_data = np.genfromtxt(self.point_file,usecols=(1,0,2,3,4))

        # load random locations
        if self.random_locations_fn is None:
            raise NotImplementedError("Please provide a file containing random"+\
                                      " point locations, seeding of random "+\
                                      "locations not implemented yet")
        self.points_random = np.load(self.random_locations_fn)
            
            
    def prepare_data(self):
        
        """
        prepare input data by the following steps:
             - filter point data to be within station_buffer of a station and 
               within/outside within_polygonfile and outside_polygonfile
             - project to eastings and northings
             - determine number of points and store
             - filter and project random points
             - initialise arrays to contain distances and point property data
        
        """
                
        # filter points for distance to station and within/outside specified
        # polygons
        self.point_data_filt, station_filt, polygon_filt = \
            filter_points(self.point_data,self.station_ll,self.station_buffer,
                          self.within_polygonfile,self.outside_polygonfile)
            
        # project to a projected coordinate system in metres
        self.point_data_filt_xy = np.zeros_like(self.point_data_filt)
        self.point_data_filt_xy[:,0], self.point_data_filt_xy[:,1] = \
            epsg_project(self.point_data_filt[:,0],self.point_data_filt[:,1],
                         4326,0,proj_str=self.projection)
            
        # get number of deposits
        self.n_points = len(self.point_data_filt)
            
        # make a property array
        self.make_point_property_array()
        
        # load and filter random points
        self.filter_project_random_points()
        
        # initialise distance array
        self.make_distance_array()
        
        
        
    def run_analysis(self):
        """
        
        Returns
        -------
        None.
        """
        
        for idx in range(self.gcz.size):
            t0 = time.time()
            values = self.values[:,:,idx]

            #compute 100 ohm-m contour and turn into a shapely MultiLineString
            mlinestr = compute_contour(self.gcx,
                                       self.gcy,
                                       values,
                                       self.target_value)
            pointlist = [Point(x,y) for x,y in self.point_data_filt_xy[:,:2]]
            
            # compute the distance to the contour for each gold deposit
            if mlinestr is None:
                self.distance_array['distance'][idx] = np.nan
            else:
                self.distance_array['distance'][idx] = \
                    compute_distance_to_target_contour(self.gcx,self.gcy,
                                                        values,
                                                        self.target_value,
                                                        mlinestr,
                                                        pointlist)
                
                # compute the distance from target contour for the random locations
                for r in range(self.n_repeats):
                    rpointlist = [Point(rxy) for rxy in self.points_random_xy[r]]
                    self.distance_array['distance_random'][idx,r] = \
                        compute_distance_to_target_contour(self.gcx,self.gcy,
                                                            values,
                                                            self.target_value,
                                                            mlinestr,
                                                            rpointlist)
            t1 = time.time()
            print("completed depth %1i of %1i"%(idx,self.gcz.size),t1-t0)
                                
        prefix = '%s_%1iohmm'%(os.path.basename(self.point_file)[:-4],
                                      self.target_value)
        
        self.distance_array_fn = os.path.join(self.savepath,
                                              self.prefix+'_distance_array.npy')
        self.point_property_array_fn = os.path.join(self.savepath,
                                                    prefix+'_contained_resource.npy')
        
        np.save(self.distance_array_fn, self.distance_array)
        np.save(self.point_property_array_fn, self.point_property_array)
        
        print("saved arrays",self.prefix)
            

        
    
        
        
    def filter_project_random_points(self):
        """
        Filter random points to be within station_buffer of a station and
        within/outside the polygons. Project to eastings and northings
        Returns
        -------
        None.
        """
        
        self.points_random,_,filt2 = filter_points(self.points_random, 
                                                   self.station_ll,
                                           self.station_buffer,
                                           self.within_polygonfile,
                                           self.outside_polygonfile,
                                           maxpoints=self.n_repeats*self.n_points)
        self.points_random_xy = np.zeros_like(self.points_random)
        self.points_random_xy[:,0],self.points_random_xy[:,1] = \
            epsg_project(self.points_random[:,0],self.points_random[:,1],
                         4326,0,proj_str=self.projection)
        self.points_random = self.points_random.reshape(self.n_repeats,
                                                        self.n_points,2)
        self.points_random_xy = self.points_random_xy.reshape(self.n_repeats,
                                                        self.n_points,2)
        
        
    def make_point_property_array(self):
        """
        Make and populate an array to contain point properties (filtered as for 
        x,y locations).
        Returns
        -------
        None.
        """
        
        arrdtype = [(prop,np.float) for prop in self.point_propertynames]
        self.point_property_array = np.array(np.zeros(len(self.point_data_filt)),
                                      dtype=arrdtype)        
        for i, pname in enumerate(self.point_propertynames):
            self.point_property_array[pname] = self.point_data_filt[:,i+2]
            
            
    def make_distance_array(self):
        """
        Make and populate an array to contain distances to contour.
        Returns
        -------
        None.
        """
        
        # initialise structured array to contain distance-to-contour data    
        self.distance_array = np.array(np.zeros(len(self.gcz)),
                                  dtype=[('depth',np.float),
                                          ('x',np.float,self.n_points),
                                          ('y',np.float,self.n_points),
                                          ('distance_random',np.float,(self.n_repeats,self.n_points)),
                                          ('distance',np.float,self.n_points)])
        self.distance_array['depth'] = self.gcz
        self.distance_array['x'] = self.point_data_filt_xy[:,0]
        self.distance_array['y'] = self.point_data_filt_xy[:,1]  
        
        
if __name__ == "__main__":
    wd = os.path.join('.','distance_to_contour/inputs')
    savepath = os.path.join('.','distance_to_contour/outputs')
    netcdf_file = os.path.join(wd,'test.nc')
    station_locations_file = os.path.join(wd,'test_sloc.txt')
    
    random_locations_fn = os.path.join(wd,'RandomPoints_100000points.npy')
    point_file = os.path.join(wd,'test_points.txt')
    target_value = 100
        
    distanceObj = DistanceToContour(netcdf_file, 
                                    station_locations_file,
                                    point_file, 
                                    target_value,
                                    random_locations_fn = random_locations_fn,
                                    savepath = savepath,
                                    point_propertynames=['Au','Cu','Age_Ga'],
                                    n_repeats=100)
    distanceObj.setup_and_run_analysis()