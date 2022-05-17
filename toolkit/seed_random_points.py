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
from toolkit.geophysical_model import GeophysicalModel
from toolkit.distance_to_contour_tools import epsg_project, filt_bbox
import numpy as np
from scipy.interpolate import interp1d
from mpl_toolkits.basemap import Basemap
from shapely.geometry import Point
import os


class RandomPoints(GeophysicalModel):
    
    def __init__(self, wd, **kwargs):
        """
        Class to generate randomly located points that are on land and within 
        station_buffer of station locations for a set of geophysical models.
        Inherits GeophysicalModel class

        Parameters
        ----------
        wd : string
            Working directory where files are located (outputs will also be 
                                                       saved here by default)
                                                    
        Optional Parameters
        -------------------       
        savepath : path to save to, defaults to wd
        station_buffer : buffer around stations to seed random points
        n_points : number of random points to seed
        projection : global equal area projection to seed points in. Points are
                     seeded in projected coordinates (equal area) to ensure 
                     that they are spatially even across the globe
        chunksize : points are seeded and have initial filtering done in chunks 
                    to make use of fast numpy array routines. The chunksize
                    argument determines how many points to put in each chunk
                    (note it is frequently larger than the number of points
                     as especially with resistivity models they don't cover
                     much area so only a small proportion get through initial
                     filtering)
        bounds : if known, longitude/latitude bounds can be provided. These
                 pertain to all surveys so all surveys must be within these
                 bounds. This will speed the process up
        savefilename_prefix : prefix for the filename to save the random points
                              array to.

        Returns
        -------
        None.

        """
        
        super().__init__(wd, **kwargs)
        
        self.savepath = wd
        self.n_points = 100
        self.projection = '+proj=eqearth'
        self.chunksize = 10000
        self.bounds = [-180,180,-90,90]
        self.savefilename_prefix = 'RandomPoints'
        
        # set attributes based on kwargs
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self,key,kwargs[key])
                
        # additional attribute to tell how often to output how many repeats 
        # have been completed
        self._screen_output_interval = int(10**np.round(np.log10(self.n_points/100)))
                
    def seed_random_points(self):
        """
        
        Seed random points

        Returns
        -------
        None.

        """
        
        # initialise a chunk of random points to assess
        self.compute_legal_bounds()
        
        # initialise new array to save points to
        self.random_locations = np.zeros((self.n_points,2))
        
        # read modem data and model
        self.read_models()
        self.compile_station_locations()
       
        while np.all(self.random_locations[-1] == 0):
            
            # initialise a chunk of random points
            self._initialise_random_points_chunk()
            
            # filter (roughly) to be within station area
            self._filter_bbox()
        
            # more refined filter to ensure within station_buffer of a station,
            # within station_buffer/2 of the convex hull around the stations, 
            # and on land
            self.filter_stations()
            
        # save array to file
        self.save_array()
        
            

    def compute_legal_bounds(self):
        """
        Create functions that define the legal bounds for the projection (i.e.
        between -180 and 180 degrees longitude and between -90 and 90 degrees
        latitude)

        Returns
        -------
        None.

        """
        
        # calculate eastern and western bounds so we don't return "illegal" values
        eblat = np.arange(-90,90.1,0.1)
        eblon = np.ones_like(eblat)*180
        ebx,eby = epsg_project(eblon,eblat,4326,0,proj_str=self.projection)
        
        wblat = np.arange(-90,90.1,0.1)
        wblon = np.ones_like(wblat)*-180
        wbx,wby = epsg_project(wblon,wblat,4326,0,proj_str=self.projection)
        
        #interpolation function for eastern & western bound
        self._ebf = interp1d(eby,ebx)
        self._wbf = interp1d(wby,wbx)
        
    
    def _initialise_random_points_chunk(self):
        """
        Initialise points within the bounds as stored within the attribute
        to ensure equal area representation across the globe they are seeded in 
        eastings and northings then filtered to be within legal latitude 
        longitude bounds

        Returns
        -------
        None.

        """
        
        lon0, lon1, lat0, lat1 = self.bounds
        
        # project to eastings/northings
        x0 = epsg_project(lon0,(lat0+lat1)/2,4326,0,proj_str=self.projection)[0]
        y0 = epsg_project((lon0+lon1)/2,lat0,4326,0,proj_str=self.projection)[1]
        x1 = epsg_project(lon1,(lat0+lat1)/2,4326,0,proj_str=self.projection)[0]
        y1 = epsg_project((lon0+lon1)/2,lat1,4326,0,proj_str=self.projection)[1]
        
        randeast = np.random.random(size=self.chunksize)*(x1 - x0) + x0
        randnorth = np.random.random(size=self.chunksize)*(y1 - y0) + y0
        
        filt = np.where(np.all([randeast <= self._ebf(randnorth),
                                randeast >= self._wbf(randnorth)],axis=0))
        
        randeast, randnorth = randeast[filt], randnorth[filt]
        
        randlon,randlat = epsg_project(randeast,randnorth,
                                       0,4326,proj_str=self.projection)
        if len(randlon.shape) > 1:
            randlon = randlon[0]
            randlat = randlat[0]

        self.random_locations_chunk = np.vstack([randlon,randlat]).T
        
        
    def _filter_bbox(self):
        
        """
        Filter random points to be within a rectangular box defined by station
        locations + station_buffer. This is more efficient than filtering to be
        within a distance of stations or within the survey convex hull so is
        done first to preselect random points before further filtering is
        carried out
    
        Returns
        -------
        None.

        """
        first = True
        for modObj in self.model_dict.values():
            
            minLon = modObj.slon.min() - self.station_buffer
            minLat = modObj.slat.min() - self.station_buffer
            maxLon = modObj.slon.max() + self.station_buffer
            maxLat = modObj.slat.max() + self.station_buffer
        
            if not first:
                filt2 = np.any([filt2,
                                filt_bbox(self.random_locations_chunk,
                                          minLon,maxLon,
                                          minLat,maxLat)], axis=0)
            else:
                filt2 = filt_bbox(self.random_locations_chunk,
                                  minLon,maxLon,
                                  minLat,maxLat)
                first = False
                
        self.random_locations_chunk = self.random_locations_chunk[filt2]
        
    
    def filter_stations(self):
        
        count = np.where(np.all(self.random_locations == 0,axis=1))[0][0]

        bm = Basemap()
        for i in range(self.random_locations_chunk.shape[0]):
            # print("assessing point %1i"%i)
            if count >= self.random_locations.shape[0]:
                break
            
            lon,lat = self.random_locations_chunk[i]

            # check if within station_buffer of a station, station_buffer/2
            # of the convex hull of the station, and on land.
            if Point(lon,lat).distance(self.station_multipolygon) < self.station_buffer/2:
                # print("within convex hull")
                if Point(lon,lat).distance(self.station_multipoint) < self.station_buffer:
                    # print("within station buffer")
                    if bm.is_land(lon,lat):
                        # print("on land")
                        self.random_locations[count] = [lon,lat]
                        count += 1
                        
                        if count % self._screen_output_interval == 0:
                           print("completed repeat",count,"of",self.n_points)


    def save_array(self):
        
        
        savefile = self.savefilename_prefix + '_%1ipoints'%self.n_points
        self.random_locations_fn = os.path.join(self.savepath,savefile)
        np.save(self.random_locations_fn,self.random_locations)
        

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Command line input.')
    parser.add_argument("--mst", action="store", help="Path to the Mineral-Stats-Toolkit installation directory.")
    rv = parser.parse_args()
    wd = os.path.join(rv.mst,'data/resistivity_models')
    savepath = '.'
    
    modem_dict = {
                   'Kirkby2020':['Modular_MPI_NLCG_005.rho',[3.789018683870512177e+05,6.158864711905587465e+06],
                                      '+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                      'Robertson2020':['Delamerian_AusLAMP_cropped.rho',[366832.99996643455, 6428315.2838392835],
                                        '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
    }
    
    rp = RandomPoints(wd, 
                      modem_details_dict=modem_dict,
                      savepath=savepath,
                      n_points=100000,
                      chunksize=10000000)
    rp.seed_random_points()