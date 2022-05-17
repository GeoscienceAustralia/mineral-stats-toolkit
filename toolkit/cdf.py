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
import matplotlib.pyplot as plt
from toolkit.distance_to_contour_tools import nearest_index, as_si

class CDF():
    def __init__(self, **kwargs):
        """
        

        Parameters
        ----------
        distance_array_fn : TYPE
            DESCRIPTION.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        self.distance_array_fn = None
        self.property_array_fn = None
        self.binsize = None
        self.cdf_array_fn = 'cdf'
        self.savepath = '.'
        
        for key in kwargs.keys():
            if hasattr(self,key):
                setattr(self,key,kwargs[key])
        
        if not os.path.exists(self.savepath):
            os.mkdir(self.savepath)             
        
        
        
    def read_distances(self):
        """
        Read files containing distances and properties

        Returns
        -------
        None.

        """
        self.distance_array = np.load(self.distance_array_fn)
        if self.property_array_fn is not None:
            self.property_array = np.load(self.property_array_fn)
            
            
    def filter_by_size(self):
        return
        
    def read_cdf(self):
        """

        Returns
        -------
        None.

        """
        self.cdf_array = np.load(self.cdf_array_fn)
        self.bins = self.cdf_array['bins']
        self.binc = self.cdf_array['binc']
        self.binsize = np.median(self.bins[1:] - self.bins[:-1])
        self.n_deposits = self.cdf_array['cdf'].max()

        
    def initialise_bins(self):
        if self.binsize is None:
            self.binsize = 2e3
        
        self.bins = np.arange(np.floor(self.distance_array['distance_random'].min()/\
                                      self.binsize)*self.binsize,
                              np.ceil(self.distance_array['distance_random'].max()/\
                                      self.binsize)*self.binsize,
                         self.binsize)
        self.binc = np.mean([self.bins[:-1],self.bins[1:]],axis=0)     
        
        
    def initialise_cdf_array(self):
        
        # initialise an array to contain cdf
        self.cdf_array = np.zeros(len(self.distance_array['distance']),
                             dtype=[('depth',np.float),
                                    ('binc',np.float,self.binc.shape),
                                    ('bins',np.float,self.bins.shape),
                                    ('cdf_random',np.float,
                                     (self.distance_array['distance_random'].shape[1],
                                      self.binc.shape[0])),
                                    ('cdf',np.float,self.binc.shape)
                                    ])
        self.cdf_array['depth'] = self.distance_array['depth']
        
    
    def compute_cdf(self):
        
        self.read_distances()
        self.initialise_bins()
        self.initialise_cdf_array()
        
        distances = self.distance_array['distance']
        distance_random = self.distance_array['distance_random']
        
        for ii in range(distances.shape[0]):
            
            hist_kwargs = dict(cumulative=True,
                               bins=self.bins, histtype='step',
                               density=False)
            rhist_kwargs = dict(cumulative=True,
                               bins=self.bins, histtype='step',
                               density=False)
  
            if not np.all(distances[ii]==0):
                self.cdf_array['cdf'][ii], self.cdf_array['bins'][ii], patches= \
                    plt.hist(distances[ii],**hist_kwargs)
    
            for i in range(distance_random.shape[1]):
                self.cdf_array['cdf_random'][ii,i], bins, patches = \
                    plt.hist(distance_random[ii,i],**rhist_kwargs)
            
            self.cdf_array['binc'][ii] = self.binc
        plt.close()
        self.n_deposits = self.cdf_array['cdf'].max()
        self.save_cdf_array()
        
    def compute_normed_cdf(self):
        """
        

        Raises
        ------
        AttributeError
            DESCRIPTION.

        Returns
        -------
        None.

        """

        if not hasattr(self,'cdf_array'):
            raise AttributeError("Please compute or load a CDF array first")
            

        self.cdf_normed = self.cdf_array['cdf']/self.n_deposits
        self.cdf_random_normed = np.mean(self.cdf_array['cdf_random'],axis=1)/\
            self.n_deposits
        self.cdf_random_std = np.std(self.cdf_array['cdf_random'],axis=1)/\
            self.n_deposits        

        self.dvalues = self.cdf_normed - self.cdf_random_normed
    
        self.dmax = np.nanmax(self.dvalues,axis=1)
        
        
        
    def compute_kolmogorov_smirnov(self):
        """
        Compute kolmogorov-smirnov probability that distance distribution could
        be pulled by chance

        Returns
        -------
        None.

        """
        self.compute_normed_cdf()
         
        self.kolmogorov_smirnov = np.exp(-2.*(self.n_deposits**2)*(self.dmax**2)/\
                                         (2*self.n_deposits))
        
            
    def plot_cdf(self, depth, normed=True, color='b', pad = 20e3):
        if not hasattr(self,'cdf_array'):
            raise AttributeError("Please compute or load a CDF array first")
            
        didx = nearest_index(depth, self.cdf_array['depth'])
            
        # get data to plot
        if normed:
            if not hasattr(self,'cdf_normed'):
                self.compute_normed_cdf()
            cdf, cdfr, cdfstd = self.cdf_normed[didx]*100, \
                self.cdf_random_normed[didx]*100, \
                self.cdf_random_std[didx]*100
            idx1 = np.where(cdfr==100)[0][0]
            y0, y1 = -5, 105
        else:
            cdf, cdfr, cdfstd = self.cdf_array['cdf'][didx],\
                np.mean(self.cdf_array['cdf_random'][didx],axis=1),\
                np.std(self.cdf_array['cdf_random'][didx],axis=1)
            y0, y1 = np.array([-0.05, 1.05])*self.n_deposits
            idx1 = np.where(cdfr==self.n_deposits)[0][0]
            
        idx0 = np.where(cdfr==0)[0][-1]
        
        #print(idx0, idx1)
        #print(didx)
        x0, x1 = self.cdf_array['binc'][didx,idx0] - pad, \
                 self.cdf_array['binc'][didx,idx1] + pad
            
        #print(self.cdf_array['binc'][didx,idx0])
        #print(self.cdf_array['binc'][didx,idx1])
        #print(x0)
        #print(x1)
        plt.plot(self.cdf_array['binc'][didx], cdf, color=color)
        plt.fill_between(self.cdf_array['binc'][didx],
                         cdfr - cdfstd, 
                         cdfr + cdfstd,
                         color=color,
                         alpha=0.5)
        plt.xlim(x0, x1)
                
        
    def plot_heatplot(self,ax=None,show_dmax=True,ks_depths=['min']):
        """
        Generate a heat plot showing d values (difference between )

        Parameters
        ----------
        show_dmax : TYPE, optional
            Whether or not to plot dmax (normalised difference cdf - cdf_random) 
            as a function of depth on top of the heat plot. The default is True.
        ks_depths : TYPE, optional
            Depths in metres at which to display the kolmogorov-smirnov (KS)
            probability that the distribution with respect to the contour 
            could be pulled by chance. The default is ['min'] which finds the
            depth at which the KS probability is minimised

        Returns
        -------
        None.

        """
        self.compute_normed_cdf()
        bins, depth = [self.cdf_array[att] for att in ['bins','depth']]
        if ax is None:
            ax = plt.subplot(111)
        depth = depth/1e3
        bins = bins/1e3
        
        ax.pcolor(bins,depth,self.dvalues,
                   vmin=0,
                   vmax=0.5,
                   cmap='hot_r')
        

                
                
        ax.set_xlim(-70,120)
        ax.set_ylim(150, 0)
        
        ax.set_aspect(2)
    

        if ks_depths:
            ax2 = plt.twiny()
            self.compute_kolmogorov_smirnov()
            for ks_depth in ks_depths:
                if ks_depth == 'min':
                    didx = np.where(self.kolmogorov_smirnov==np.amin(self.kolmogorov_smirnov))[0][0]
                else:
                    didx = nearest_index(ks_depth,depth)
                pkslabel = r"P$_{KS}=$"+"\n${0:s}$".format(as_si(self.kolmogorov_smirnov[didx]))
                ax2.text(self.dmax[didx],depth[didx],pkslabel,
                          ha='left',va='center',color='k',fontsize=10)            

        if show_dmax:
            if not ks_depths:
                ax2 = plt.twiny()
            ax2.plot(self.dmax,depth,'k-')
            ax2.set_xlim(0,np.around(self.dmax.max()+0.2,1))
            ax2.set_ylim(150, 0)       
        
    def show(self):
        plt.show()
        
        
        
    def save_cdf_array(self):
        
        np.save(os.path.join(self.savepath,self.cdf_array_fn),
                self.cdf_array)
        
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Command line input.')
    parser.add_argument("--mst", action="store", help="Path to the Mineral-Stats-Toolkit installation directory.")
    rv = parser.parse_args()
    wd = os.path.join(rv.mst,'data/outputs')
    savepath = '.'
    
    cdf_obj = CDF(distance_array_fn=os.path.join(wd,'test_points_100ohmm_distance_array.npy'),
                  cdf_array_fn = 'cdf',
                  savepath=savepath,
                  binsize=2e3
                  )
    
    cdf_obj.compute_cdf()

    cdf_obj.plot_cdf(50e3)
                     