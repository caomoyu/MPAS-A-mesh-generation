 
import numpy as np
import argparse
import os
import xarray
import jigsawpy as jig
import jigsaw_util as jutil
from mpas_tools.mesh.creation.jigsaw_to_netcdf import jigsaw_to_netcdf
from mpas_tools.mesh.conversion import convert
from mpas_tools.io import write_netcdf
from mpas_tools.ocean.inject_meshDensity import inject_spherical_meshDensity
earth_radius=6371.0e3
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point
import os
import shutil


res1 = 0.5
res2 = 1
res3 = 10
resall = 30

# clon, clat = 114,22
clon, clat = 114.17,22.29
# dlat = 0.5 / 30  # Smaller lat-lon grid interval for finer resolution near equator
# dlon = dlat
dlat=0.5/30
dlon=0.5/30
nlat = int(180.0 / dlat) + 1
nlon = int(360.0 / dlon) + 1

lat = np.linspace(-90., 90., nlat)
lon = np.linspace(-180., 180., nlon)
lons, lats = np.meshgrid(lon, lat)

# Calculate distances to center (lat=clat, lon=clon)
earth_radius_km = 6371.0
clat_rad, clon_rad = np.radians(clat), np.radians(clon)
lat_rad, lon_rad = np.radians(lats), np.radians(lons)

delta_sigma = np.arccos(
    np.sin(clat_rad) * np.sin(lat_rad) +
    np.cos(clat_rad) * np.cos(lat_rad) * np.cos(lon_rad - clon_rad)
)
dists = earth_radius_km * delta_sigma  # Distance in kilometers

# Initialize cell width with the largest (100 km) resolution
cellWidth = np.full_like(dists, 30)
# Apply 4 km resolution within the inner 100 km radius
# cellWidth[dists <= 30] = 0.5

# Apply 10 km resolution within 600 km, transitioning from 4 km to 10 km
# transition_zone_10km = (dists > 30) & (dists <= 35)
# slope_4_to_10 = (1 - 0.5) / (35 - 30)
# cellWidth[transition_zone_10km] = 0.5 + slope_4_to_10 * (dists[transition_zone_10km] - 30)

cellWidth[dists <= 50] = 0.5
cellWidth[(dists > 50) & (dists <= 200)] = 1
# cellWidth[(dists > 200) & (dists <= 600)] = 10
# mask = (lons >= 95) & (lons <= 140) & (lats >= -10) & (lats <= 25)
# cellWidth[~mask] = 120  # Set a high cell width outside the box to effectively mask them


# sea_polygon = Polygon([(-10, 95), (25, 95), (25, 140), (-10, 140)])
# for i in range(lat.size):
#     for j in range(lon.size):
#         point = Point(lon[j], lat[i])
#         if sea_polygon.contains(point):
#             cellWidth[i, j] = 10

transition_zone_10km = (dists > 200) & (dists <= 1000)
slope_4_to_10 = (10 - 1) / 800
cellWidth[transition_zone_10km] = 1 + slope_4_to_10 * (dists[transition_zone_10km] - 200)

# mask = (lons >= 90) & (lons <= 150) & (lats >= -20) & (lats <= 30)
# cellWidth[~mask] = 120  # Set a high cell width outside the box to effectively mask them


# cellWidth[(dists > 1000) & (dists <= 600)] = 20
transition_zone_10km = (dists > 1000) & (dists <= 1200)
slope_4_to_10 = (30-10) / (200)
cellWidth[transition_zone_10km] = 10 + slope_4_to_10 * (dists[transition_zone_10km] - 1000)


out_dir='./grid_test24'
out_base='./grid_test24'
out_basepath=out_dir+"/"+out_base
out_filename=out_dir+"/"+out_base+".mpas.nc"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)
else:
    print("Base dir already exists: ", out_dir)


mesh_file = jutil.jigsaw_gen_sph_grid(cellWidth, lon, lat, basename=out_basepath) 


jigsaw_to_netcdf(msh_filename=mesh_file,
                        output_name=out_basepath+'_triangles.nc', on_sphere=True,
                        sphere_radius=1.0)

#convert to mpas grid specific format
write_netcdf(
        convert(xarray.open_dataset(out_basepath+'_triangles.nc'), 
        dir=out_dir,
        graphInfoFileName=out_basepath+"_graph.info"),
        out_filename)


inject_spherical_meshDensity(cellWidth, lon, lat,
                                     mesh_filename=out_filename)

shutil.copy('/disk/r106/ycaobx/MPAS-Tools/MPAS-BR/mesh_cye_SEA2.py', \
    out_dir+'/'+'mesh_cye_SEA2.py')