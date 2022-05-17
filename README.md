# mineral-stats-toolkit

This repository contains the code to analyse the distance from point data to a contour in a geoscientific dataset, compared to a set of random locations, as described in Kirkby et al (in Review). Users of this repository are asked to cite Kirkby et al (in Review).

Kirkby, A., Czarnota, K., Huston, D., Champion, D., Doublier, M., Bedrosian, P., Duan, J., Heinson, G., 2022. Lithospheric conductors reveal source regions of convergent margin mineral systems. Scientific Reports, 12, 8190, doi: https://doi.org/10.1038/s41598-022-11921-2

## Dependencies

* mtpy
* netCDF4
* geopandas
* mpl_toolkits.basemap - this is depreciated and should be replaced with cartopy

## Running

The repository can be installed simply by cloning from Github:

```
git clone https://github.com/GeoscienceAustralia/mineral-stats-toolkit.git
```

After adding the installation directory to your PYTHONPATH environment variable, you can verify the codes by running the included examples. For instance:

```
cd ..
mkdir tmp
cd tmp
python ../mineral-stats-toolkit/examples/run.py --mst ../mineral-stats-toolkit
python ../mineral-potential-toolkit/examples/plot_cdf.py
```

## Input data

The user needs to specify two input values:
- target value (e.g. resistivity value) of interest to contour
- number of repeats to run when computing the random locations

The following input files are also required (see data folder):
- a netCDF file containing a 3D numpy array containing values (dimensions m * n * o) and 1d arrays containing latitude (m), longitude (n) and depth (o). No data values denoted by nans in the numpy array.
- a column-based text file with station locations
- a point dataset (e.g. mineral deposit locations)
- a .npy file (standard numpy binary format) containing random locations covering the same spatial area as the data in the netcdf file

There are functions in the repository to convert a collection of ModEM format resistivity models to the netCDF format required, and to generate the random points.

## Output data

The analysis generates three outputs:
- File containing the x, y, locations of each point and its distance to the contour
- File containing the CDFs for the true points and for random locations.

## Algorithms

The key algorithms are demonstrated in the examples folder:
create_inputs_from_modem.py generates a netCDF file from a collection of resistivity models, station locations file and .npy file with random locations
run.py runs the analysis:
- first, it computes the distance to contour
- then, it computes CDFs for the random and true point locations.
Visualisation tools are under development.

