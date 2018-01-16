import numpy as np
import matplotlib.pyplot as plt
import cmocean
from mpl_toolkits.basemap import Basemap
import netCDF4
from scipy import interpolate
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.patheffects as PathEffects
from matplotlib.path import Path
import urllib2
import warnings
import numpy.ma as ma
import time

warnings.filterwarnings("ignore",".*invalid value encountered in less.*")

version_number = 1.1
'''
Version 1.0 - first working draft
1.1 - added analysis time global attribute, attempted to make upper levels
randomly decrease (doesn't work yet)
'''

my_nextcloud = '/Users/briangreene/Nextcloud/thermo/'

# Oklahoma grid
llcrnrlat = 33.5
urcrnrlat = 37.2
llcrnrlon = -103.2
urcrnrlon = -94.0

# Grid spacing - degrees lat, lon
gridspace = 0.01

# Initialize map
fig_map, ax_map = plt.subplots(1, figsize=(12, 6.75))
m = Basemap(projection='merc', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, 
    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, resolution='l', ax=ax_map)

# Draw patches over not Oklahoma
shapefile = my_nextcloud + 'documentation/States/states_21basic/states'
m.readshapefile(shapefile, 'states', drawbounds = False)
patch = []
patch_ok = []
for info, shape in zip(m.states_info, m.states):
	if info['STATE_NAME'] != 'Oklahoma':
		patch.append( Polygon(np.array(shape), True) )
	else:
		patch_ok.append( Polygon(np.array(shape), True))

# get OK border polygon points
xy = patch_ok[0].get_xy()

# dimensions of new meshgrid
xnew = np.arange(llcrnrlon, urcrnrlon, gridspace)
ynew = np.arange(llcrnrlat, urcrnrlat, gridspace)
xx, yy = np.meshgrid(xnew, ynew)
# flatten
xxyyflat = np.dstack((xx, yy))
xxyyflat = xxyyflat.reshape(-1, 2)

URL = 'http://www.mesonet.org/data/public/mesonet/current/current.csv.txt'
fd = urllib2.urlopen(URL)
data_long = fd.read()
data = data_long.split('\n')
lat_, lon_, yr_, mo_, da_, hr_, mi_, T_, Temps, lats, lons, yr, mo, da, hr, mi \
	= ([] for i in range(16))
for i in np.arange(1, len(data)-1):
	data_short = data[i].split(',')
	lat_.append(float(data_short[3]))
	lon_.append(float(data_short[4]))
	yr_.append(float(data_short[5]))
	mo_.append(float(data_short[6]))
	da_.append(float(data_short[7]))
	hr_.append(float(data_short[8]))
	mi_.append(float(data_short[9]))
	T_.append(data_short[10])

for i in np.arange(len(T_)):
	if T_[i] != ' ':
		Temps.append(float(T_[i]))
		lats.append(lat_[i])
		lons.append(lon_[i])
		yr.append(yr_[i]) 
		mo.append(mo_[i]) 
		da.append(da_[i]) 
		hr.append(hr_[i]) 
		mi.append(mi_[i]) 

timeValid_str = '%04d%02d%02d_%02d%02d' % (yr[0], mo[0], da[0], hr[0], mi[0])

# transopose to basemap coordinates
lon, lat, xplot, yplot = m.makegrid(len(xnew), len(ynew), returnxy=True)
x1, y1 = m(lons, lats)

# Extend values to perimiter
def nearest(Lon, Lat):
	# returns index of closest mesonet site
	dist = np.sqrt((Lon - np.array(x1))**2. + (Lat - np.array(y1))**2.)
	return np.argmin(dist)

Tperim, latperim, lonperim = ([] for i in range(3))
for i in range(len(xy)):
	iclose = nearest(xy[i, 0], xy[i, 1])
	Tperim.append(Temps[iclose])
	latperim.append(xy[i, 1])
	lonperim.append(xy[i, 0])

x2 = x1 + lonperim
y2 = y1 + latperim
T2 = Temps + Tperim

# mesonet lat lon pairs in basemap coords
points = np.array(zip(x2, y2))

# interpolate
t0 = time.time()
Tnew = interpolate.griddata(points, T2, (xplot, yplot), method='cubic')
print "Interpolate time: %s seconds" % (time.time() - t0)
# Determine which points inside OK
mpath = Path(xy)
xxyyplot = np.dstack((xplot, yplot))
xxyyplot_flat = xxyyplot.reshape(-1, 2)
tstart = time.time()
# use invert because this shapefile goes clockwise, mpath assumes ccw
mask_flat = np.invert(mpath.contains_points(xxyyplot_flat))
mask = mask_flat.reshape(xplot.shape)
print "Mask time: %s seconds" % (time.time() - tstart)

# Cast mask as nans - this is sfc level
T0 = Tnew
T0[mask] = np.nan

# Create a few new levels
T1 = T0 - (5. * np.random.random_sample(T0.shape))
T2 = T1 - (5. * np.random.random_sample(T0.shape))
T3 = T2 - (5. * np.random.random_sample(T0.shape))

'''
Now we have interplated data on specified grid size and spacing, with a few
levels chosen. Begin creating nc file.
'''

# Initialize nc file
fname = '/Users/briangreene/Desktop/sample3Dmeso.nc'
rootgrp = netCDF4.Dataset(fname, 'w', format='NETCDF4')
# Initialize dimensions - level, lat, lon
level = rootgrp.createDimension('level', None)
lat = rootgrp.createDimension('lat', len(ynew))
lon = rootgrp.createDimension('lon', len(xnew))
x = rootgrp.createDimension('x', None)
# Initialize variables
# Dimensional info
latitude = rootgrp.createVariable('lat', 'f4', ('lat',))
longitude = rootgrp.createVariable('lon', 'f4', ('lon',))
levels = rootgrp.createVariable('level', 'i4', ('level',))
x_station = rootgrp.createVariable('x', 'f4', ('x',))
y_station = rootgrp.createVariable('y', 'f4', ('x',))
# Data
T = rootgrp.createVariable('Temperature', 'f4', ('lat','lon','level'))
Tstation = rootgrp.createVariable('Station_Temp', 'f4', ('x',))
# Global atts
rootgrp.description = 'Sample 3D Mesonet nc file'
rootgrp.version = version_number
rootgrp.lllat = llcrnrlat
rootgrp.urlat = urcrnrlat
rootgrp.lllon = llcrnrlon
rootgrp.urlon = urcrnrlon
rootgrp.timeValid = timeValid_str
# Variable attributes
T.units = 'C'
T.name_long = 'Temperature'
# Assign Values
latitude[:] = ynew
longitude[:] = xnew
T[:, :, 0] = T0
T[:, :, 1] = T1
T[:, :, 2] = T2
T[:, :, 3] = T3
Tstation[:len(Temps)] = Temps
x_station[:len(Temps)] = lons
y_station[:len(Temps)] = lats
# Assign to unlimited dimension variable
levels = range(4)
x = range(len(Temps))

# Close file
rootgrp.close()

print '>>File created successfully'
