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
import numpy as np
from toolkit.distance_to_contour_tools import epsg_project
from netCDF4 import Dataset
import os
from scipy.interpolate import RegularGridInterpolator




class MergeModels(GeophysicalModel):
    
    
    def __init__(self, wd, **kwargs):
        """
        Class to manage combining geophysical models into a single
        grid and writing to a netcdf file. 
        
        Models will be combined in the order provided in the input dict
        with no smoothing or interpolation applied between edges if there is 
        overlap between models (i.e. resistivity values from subsequent models 
        will overwrite previous values)

        Parameters
        ----------
        wd : str
            Working directory where subfolders with ModEM inversion files are 
            located.

                
        Optional Parameters
        -------------------
        savepath : path to save to, defaults to wd
        cellsize : cellsize in netcdf file to save values at
        depth_list : depth values to write resistivity values to in netcdf
        stationxy_fn_basename : name of station xy files (same for all 
            inversions). Code will search for station xy file within the 
            subfolder defined by modem_details_dict. Contains 2 columns, 
            longitude and latitude for station xy points
        filename_prefix : prefix for saving. Files will be called prefix.nc
            (netcdf file) and prefix_sloc.txt (station locations).




        
        Attributes
        ----------
        

        Returns
        -------
        None.
        
        
        Example
        -------
        

        """
        super().__init__(wd, **kwargs)
        
        self.wd = wd
        
        # set defaults
        self.savepath = wd
        self.cellsize = 0.1
        self.depth_list = np.arange(0,151e3,2e3)
        

        
        
        self.savefilename_prefix = 'ResModel'
        
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self,key,kwargs[key])
        
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)        
        
        
    def merge_models(self):
        """
        Read in ModEM models and interpolate onto a single grid.

        Returns
        -------
        None.

        """
        self.read_models()
        self.compile_station_locations()
        self.interpolate_to_grid()

        
    def create_netcdf(self):
        """
        Read in ModEM models, interpolate to grid, and save to netcdf file.
        Save stations as a separate textfile.

        Returns
        -------
        None.

        """
        
        self.merge_models()
        self.write_netcdf()
        
    
      
    def interpolate_to_grid(self):
        """
        Interpolate 

        Returns
        -------
        None.

        """
        
        # make a master grid to store data in
        minLon,minLat,maxLon,maxLat = -180,-90,180,90
        
        self.longitude = np.arange(minLon,maxLon+1e-8,self.cellsize)
        self.latitude = np.arange(minLat,maxLat+1e-8,self.cellsize)
        self.values = np.zeros((self.longitude.size,
                                     self.latitude.size,
                                     self.depth_list.size))
        
        # dictionary to store min/max lat/lon in each model
        self.bounds_dict = {}
        
        
        for key in self.model_dict.keys():
            modObj = self.model_dict[key]
            proj_str = self.modem_details_dict[key][2]
        
            # get min/max longitude and latitude
            minLon,maxLon,minLat,maxLat = [modObj.slon.min(),
                                           modObj.slon.max(),
                                           modObj.slat.min(),
                                           modObj.slat.max()]
        
            # force onto even lat/lon intervals
            minLon,minLat = [np.floor(val-self.station_buffer) \
                             for val in [minLon, minLat]]
            maxLon,maxLat = [np.ceil(val+self.station_buffer) \
                             for val in [maxLon, maxLat]]

            # define longitude and latitude to interpolate to
            loni = np.arange(minLon,maxLon+self.cellsize,self.cellsize)
            lati = np.arange(minLat,maxLat+self.cellsize,self.cellsize)
            loni,lati,zi = np.meshgrid(loni,lati,self.depth_list)            
            
            print('projecting locations....')
            
            xi,yi = epsg_project(loni,lati,
                                 4326,0,proj_str=proj_str)

            print('creating interpolation function....')
            resfunc = RegularGridInterpolator((modObj.gcy,
                                               modObj.gcx,
                                               modObj.gcz),
                                               np.log10(modObj.res_model),
                                               bounds_error=False)

            
            print('creating points to interpolate to...')
            
            xyzi = np.vstack([yi.flatten(),xi.flatten(),zi.flatten()]).T
            ny,nx,nz = xi.shape

            print('interpolating to grid....')
            

            resi = 10**resfunc((xyzi)).reshape(ny,nx,nz).transpose(1,0,2)
        
            
            # mesh into master grid
            
            i0 = np.where(np.abs(loni[0,0,0]-self.longitude) < 1e-8)[0][0]
            i1 = i0 + loni.shape[1]
            
            # print(lati,self.latitude)
            j0 = np.where(np.abs(lati[0,0,0]-self.latitude) < 1e-8)[0][0]
            j1 = j0 + lati.shape[0]    
            
            self.values[i0:i1,j0:j1][np.isfinite(resi)] = resi[np.isfinite(resi)]
            
        # trim to the extent of the data + 5 padding cells on each side
        ivals = np.where(self.values.sum(axis=1))[0]
        jvals = np.where(self.values.sum(axis=0))[0]
        i0,i1 = ivals[0]-5,ivals[-1]+5
        j0,j1 = jvals[0]-5,jvals[-1]+5
        
        self.longitude = self.longitude[i0:i1]
        self.latitude = self.latitude[j0:j1]
        self.values = self.values[i0:i1,j0:j1]
        self.values[self.values==0] = np.nan    
            
        
        self.values = self.values.transpose(1,0,2)
            
        
    def write_netcdf(self):
        """
        Save resistivity model to a netcdf file and station locations to a
        text file.

        Returns
        -------
        None.

        """
        self.netcdf_fn = os.path.join(self.savepath,
                                  self.savefilename_prefix+'.nc')

        with Dataset(self.netcdf_fn,
                     'w',format='NETCDF4') as ds:
            ds.createDimension('latitude',self.latitude.size)
            ds.createDimension('longitude',self.longitude.size)
            ds.createDimension('depth',self.depth_list.size)
            
            
            latitude = ds.createVariable('latitude',float,('latitude'))
            longitude = ds.createVariable('longitude',float,('longitude'))
            
            depth = ds.createVariable('depth',float,('depth'))
            resistivity = ds.createVariable(self.propertyname,float,
                                            ('latitude','longitude','depth'))
        
            
            
            latitude[:] = self.latitude
            longitude[:] = self.longitude
            
            depth[:] = self.depth_list
            resistivity[:,:,:] = self.values
            
        self.station_locations_fn = os.path.join(self.savepath,
                                        self.savefilename_prefix+'_sloc.txt')
        np.savetxt(self.station_locations_fn,
                   self.station_ll)       
        
        

        
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Command line input.')
    parser.add_argument("--mst", action="store", help="Path to the Mineral-Stats-Toolkit installation directory.")
    rv = parser.parse_args()
    wd = os.path.join(rv.mst,'data/resistivity_models')
    savepath = '.'

    depth_list = np.hstack([np.arange(0,60,2), 
                            np.arange(60,150,5)])*1e3
    
    modem_details_dict = {
                    'Robertson2020':['Delamerian_AusLAMP_cropped.rho',[366832.99996643455, 6428315.2838392835],
                                    '+proj=utm +zone=54 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                    'Kirkby2020':['Modular_MPI_NLCG_005.rho',[3.789018683870512177e+05,6.158864711905587465e+06],
                                         '+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '],
                   }
    
    
    resmodel = MergeModels(wd,
                           savepath = savepath,
                           station_buffer=0.7,
                           modem_details_dict=modem_details_dict,
                           savefilename_prefix='test',
                           depth_list=depth_list)
    resmodel.create_netcdf()