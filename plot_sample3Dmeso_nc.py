import time
t0 = time.time()
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
import cmocean
from mpl_toolkits.basemap import Basemap
import netCDF4
import warnings
import numpy.ma as ma
print 'Import time: %f seconds' % (time.time() - t0)

'''
Created 27 December 2017
Plots sample 3D mesonet nc file
Select level in user parameters
'''

###################
# User Parameters #
###################

plotlevel = 0 # lowest = 0

###############
# Import data #
###############

fname = '/Users/briangreene/Desktop/sample3Dmeso.nc'
df = netCDF4.Dataset(fname)

T = df.variables['Temperature'][:, :, plotlevel]
lat = df.variables['lat'][:]
lon = df.variables['lon'][:]
station_lons = df.variables['x'][:]
station_lats = df.variables['y'][:]
Tstation = df.variables['Station_Temp'][:]
timeValid = df.timeValid

# Mask nans
Tplot = ma.masked_where(np.isnan(T), T)


##################
# Initialize map #
##################

# Corners
llcrnrlat = df.lllat
urcrnrlat = df.urlat
llcrnrlon = df.lllon
urcrnrlon = df.urlon

# Figure properties
t1 = time.time()
fig_map, ax_map = plt.subplots(1, figsize=(12, 6.75))
m = Basemap(projection='merc', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, 
    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, resolution='i', ax=ax_map)
print 'Time creating basemap: %s seconds' % (time.time() - t1)
# Meshgrid for plotting
a, b, xplot, yplot = m.makegrid(len(lon), len(lat), returnxy=True)
x1, y1 = m(station_lons, station_lats)

####################
# Convert Datetime #
####################
dtValid = datetime.strptime(timeValid, '%Y%m%d_%H%M')
dtValid_str = dtValid.strftime('%d %b %Y %H:%M')

########
# Plot #
########

# Contour fill
# vmin = np.nanmin(T)
# vmax = np.nanmax(T)
vmin = -10.
vmax = 110.

cfax = m.pcolormesh(xplot, yplot, Tplot, cmap=plt.get_cmap('nipy_spectral'), 
	vmin=vmin, vmax=vmax)
# cfax = m.pcolormesh(xplot, yplot, Tplot, cmap=cmocean.cm.thermal, 
# 	vmin=vmin, vmax=vmax)
cbar = m.colorbar(cfax)

# Freezing line
Tfreeze = np.array([32.0])
cax1 = m.contour(xplot, yplot, Tplot, levels=Tfreeze, colors='r')

# Station values
for i in np.arange(len(Tstation)):
	plotstring = '%d' % Tstation[i]
	txt= plt.text(x1[i], y1[i], plotstring, fontweight='bold', 
		ha='center', va='center')
	txt.set_path_effects([PathEffects.withStroke(linewidth=3, 
		foreground='w')])

# Labels
title_str = 'Oklahoma Mesonet Temerature - %s CST - Level %d' % (dtValid_str, 
	plotlevel)
ax_map.set_title(title_str)
cbar.ax.set_ylabel('Temperature ($^\circ$F)')

m.drawcounties()
m.drawstates()
print 'Total run time: %s seconds' % (time.time() - t0)
plt.show()

