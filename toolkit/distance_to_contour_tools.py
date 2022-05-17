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

from shapely.geometry import LineString, MultiLineString, shape, Point
import matplotlib.pyplot as plt
import fiona
import numpy as np
from scipy.spatial.distance import cdist
from scipy.interpolate import RegularGridInterpolator, NearestNDInterpolator

def compute_contour(gcx_grid,gcy_grid,values,target_value, pad=1,
                    return_linestr_list=False):
    """
    create a shapely MultiLineString by contouring 2D gridded data.

    Parameters
    ----------
    gcx_grid : numpy array
        2d array containing x coordinates of parameter values.
    gcy_grid : numpy array
        2d array containing y coordinates of parameter values
    values : numpy array
        2d array containing values to contour.
    target_value : float or int
        contour value in values array
    pad : int, optional
        number of cells to clip off each side of values, gcx_grid, and gcy_grid
        arrays. The default is 0.
    return_linestr_list : bool, optional
        return the contour as a list of shapely linestrings as well as a 
        multilinestring. The default is False.

    Returns
    -------
    shapely multilinestring
        contour set as a shapely multilinestring.

    """
    plt.figure()
    cs = plt.contour(gcx_grid[pad:-pad,pad:-pad],
                     gcy_grid[pad:-pad,pad:-pad],
                     values[pad:-pad,pad:-pad],levels=[target_value])
    plines = cs.collections[0].get_paths()
    linestr_list = []
    for pline in plines:
        vert = pline.vertices
        x = vert[:,0]
        y = vert[:,1]
        linestr_xy = [(ii[0], ii[1]) for ii in zip(x,y)]
        linestr_list.append(linestr_xy)

    if len(linestr_list) > 0:
        if len(linestr_list[0]) > 1:
            mlinestr = MultiLineString(linestr_list)
        else:
            mlinestr = None
    else:
        mlinestr = None
    plt.close()
    
    if return_linestr_list:
        return mlinestr,linestr_list
    else:
        return mlinestr
    
    
def epsg_project(x, y, epsg_from, epsg_to, proj_str=None):
    """
    project some xy points using the pyproj modules

    Parameters
    ----------
    x : integer or float
        x coordinate of point
    y : integer or float
        y coordinate of point
    epsg_from : int
        epsg code of x, y points provided. To provide custom projection, set
        to 0 and provide proj_str
    epsg_to : TYPE
        epsg code to project to. To provide custom projection, set
        to 0 and provide proj_str
    proj_str : str
        Proj4 string to provide to pyproj if using custom projection. This proj
        string will be applied if epsg_from or epsg_to == 0. 
        The default is None.

    Returns
    -------
    xp, yp
        x and y coordinates of projected point.

    """
    
    try:
        import pyproj
    except ImportError:
        print("please install pyproj")
        return
    
    
    proj_list = []
    for epsg in [epsg_from,epsg_to]:
        if epsg == 0:
            proj_list.append(pyproj.Proj(proj_str))
        else:
            proj_list.append(pyproj.Proj(init='epsg:%1i'%epsg))

    p1, p2 = proj_list
    

    return pyproj.transform(p1, p2, x, y)



def shapefile_list_to_polygon_dict(outside_polygonfile,
                                   within_polygonfile):
    """
    
    turn a list of shapefiles into a polygon dict with keys outside_polygonfile,
    within_polygonfile

    Parameters
    ----------
    outside_polygonfile : str
        outside_polygon shapefile .
    within_polygonfile : str
        within_polygon shapefile.
        
    Returns
    -------
    polygon_dict : dict
        dictionary containing list of polygons under keys outside_polygonfile 
        and within_polygonfile.
        
    """
    
    
    # read murray basin shapefile
    polygon_dict_keys = ['outside_polygonfile','within_polygonfile']
    polygon_dict = {}

    for i,polygonfile in enumerate([outside_polygonfile,within_polygonfile]):
        if polygonfile is not None:
            key = polygon_dict_keys[i]
            
            polygon_dict[key] = []

            sf = fiona.open(polygonfile)

            for feature in sf:
                polygon_dict[key].append(shape(feature['geometry']))

                        
    return polygon_dict


def check_within_poly(rxy,polygon_list,threshold_list):
    
    within_poly = False
    for ii,polygon in enumerate(polygon_list):
        if polygon.distance(Point(rxy)) < threshold_list[ii]:
            within_poly = True
            break
        elif Point(rxy).within(polygon):
            within_poly = True
            break
        

def check_inside_outside_polygon(rxy,polygon_dict):
    """
    check if a point is inside inside_polygon and outside outside_polygon
    with distance threshold given by threshold_dict

    Parameters
    ----------
    rxy : list or array
        x,y coordinates of a single point.
    polygon_dict : dict
        dictionary containing outside_polygonfile and within_polygonfile as 
        keys
    threshold_dict : dict
        distance thresholds to apply, the point needs to be this number of 
        metres inside or outside the polygon

    Returns
    -------
    append : bool
        True/False value as to whether the point is inside inside_polygon and
        outside outside_polygon with the given points.

    """

    
    
    
    if len(polygon_dict) == 0:
        return True
    else:
        append = True
        if 'within_polygonfile' in polygon_dict.keys():
            # first, check if it is outside the outside_polygon
            # check it is within the 'within_poly'
            append = check_within_poly(rxy,polygon_dict['within_polygonfile'])
    
        if 'outside_polygonfile' in polygon_dict.keys():
            if append:
                append = not check_within_poly(rxy,polygon_dict['outside_polygonfile'])
                
        return append


def filter_points(points,station_xy,station_buffer,within_polygonfile,
                  outside_polygonfile, maxpoints=None):
    """
    
    filter points to be within a certain buffer of a station and within/outside
    polygons as defined by polygon_dict

    Parameters
    ----------
    point_xy : numpy array, shape (n, 2)
        x, y locations of points to filter.
    station_xy : numpy array, shape (m, 2)
        x, y locations of stations
    station_buffer : float
        buffer around stations to retain points
    within_polygonfile : str
        shapefile with polygon/s covering the areas in which to include points.
        points outside the polygon are exluded.
    outside_polygonfile : str
        shapefile with polygon/s covering the areas in which to exclude points.
        points inside the polygon are exluded.    

    Returns
    -------
    None.

    """
    
    dist_to_station = np.amin(cdist(station_xy,points[:,:2]),axis=0)
    station_filt = dist_to_station < station_buffer
    points = points[station_filt]
    
    polygon_dict = shapefile_list_to_polygon_dict(outside_polygonfile,
                                                  within_polygonfile)
    # make a new pointlist and only append if within buffer of station 
    # and meets criteria in terms of polygons
    polygon_filt = np.zeros(len(points),dtype=bool)

    count = 0
    for i,xy in enumerate(points[:,:2]):
        if check_inside_outside_polygon(xy,polygon_dict):
            polygon_filt[i] = True
            count += 1
            if maxpoints is not None:
                if count >= maxpoints:
                    break
            
    points = points[polygon_filt]
    
    return points, station_filt, polygon_filt


def compute_distance_to_target_contour(gcx,gcy,values,target_value,
                                       contour_mlinestr,pointlist):

    # resistivity interpolation function in log space
    if len(gcx.shape) == 1:
        fres = RegularGridInterpolator((gcx,gcy),
                                       np.log10(values.T),
                                       method='linear')
    else:
        gcx = gcx[np.isfinite(values)]
        gcy = gcy[np.isfinite(values)]
        values = values[np.isfinite(values)]
        
        fres = NearestNDInterpolator(np.vstack([gcx.flatten(),
                                                gcy.flatten()]).T,
                                     np.log10(values).flatten())
        
    # get xy locations of filtered stations
    xy_list = np.array([[pp.coords.xy[0][0],
                         pp.coords.xy[1][0]] for pp in pointlist])

    # compute distance
    distances = np.array([pp.distance(contour_mlinestr) for pp in pointlist])
    

    # get resistivity value at station
    resval_list = 10.**np.array(fres(xy_list))

    
    # set distance for resistivities < target (i.e. within contour) to -distance
    distances[resval_list < target_value] = -distances[resval_list < target_value]

    
    return distances
    
def filt_bbox(arr,minLon,maxLon,minLat,maxLat):
    return np.all([arr[:,0] <= maxLon,
                   arr[:,0] >= minLon,
                   arr[:,1] <= maxLat,
                   arr[:,1] >= minLat],axis=0)


def nearest_index(val,array):
    """
    find the index of the nearest value in the array
    :param val: the value to search for
    :param array: the array to search in
    
    :return: index: integer describing position of nearest value in array
    
    """
    # absolute difference between value and array
    diff = np.abs(array-val)
    
    return np.where(diff==min(diff))[0][0]


def as_si(x, ndp=0):
    """
    format a number in scientific notation i.e. x 10 -3 etc
    modified from https://stackoverflow.com/questions/31453422/displaying-numbers-with-x-instead-of-e-scientific-notation-in-matplotlib/31453961

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    ndp : TYPE, optional
        DESCRIPTION. The default is 0.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times10^{{{e:d}}}'.format(m=m, e=int(e))