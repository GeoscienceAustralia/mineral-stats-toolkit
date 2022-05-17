# mineral-stats-toolkit

This repository contains the code to analyse the distance from point data to a contour in a geoscientific dataset, compared to a set of random locations, as described in Kirkby et al (2022). Users of this repository are asked to cite Kirkby et al (2022).

Kirkby, A., Czarnota, K., Huston, D., Champion, D., Doublier, M., Bedrosian, P., Duan, J., Heinson, G., 2022. Lithospheric conductors reveal source regions of convergent margin mineral systems. Scientific Reports, 12, 8190, doi: https://doi.org/10.1038/s41598-022-11921-2

## Dependencies

* mtpy
* netCDF4
* geopandas
* mpl_toolkits.basemap - this is depreciated and should be replaced with cartopy

## Running

The example .rho files in the scripts can be downloaded from:
http://dx.doi.org/10.26186/131889 (Kirkby2020)
http://pid.geoscience.gov.au/dataset/ga/144077 (Robertson2020)

The example deposit dataset can be downloaded from the supplementary data of Kirkby et al. (2022).

In order to run the analysis, make some directories to contain inputs and outputs (e.g. C:/tmp/inputs and C:/tmp/outputs). Place the resistivity models (.rho files) and deposit dataset (the example uses Orogenic_gold.txt) in C:/tmp/inputs.

The scripts should be run in the following order:
#### 1. create_inputs_from_modem.py
Creates a netcdf file with merged resistivity models (resmodels.nc), a station locations file (resmodels_sloc.txt) and randomly located points covering the area (within 0.7 degrees of latitude/longitude of a station) covered by the resistivity models
#### 2. run.py
Creates .npy files containing the distance to the specified contour ('orogenic_gold_100ohmm_distance_array.npy') and size of each deposit ('Orogenic_gold_100ohmm_contained_resource.npy')
#### 3. compute_cdf.py
Computes a cumulative distribution function for the distances to contour
#### 4. plot_cdf.py
Plots a cumulative distribution function for deposits and random, and a "heat plot" as shown in Figure 2 of Kirkby et al. 2022 which shows the difference, D, between deposits and random (on the cdf plot) as a colour image as a function of depth and distance to contour.


## Output data

The analysis generates three outputs:
- File containing the x, y, locations of each point and its distance to the contour
- File containing the CDFs for the true points and for random locations.

## Algorithms

The key algorithms are demonstrated in the examples folder:
create_inputs_from_modem.py generates a netCDF file from a collection of resistivity models, station locations file and .npy file with random locations
run.py runs the analysis, by computing the distance to contour. compute_cdf.py then computes CDFs for the random and true point locations.

