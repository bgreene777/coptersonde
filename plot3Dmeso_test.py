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

# Oklahoma
llcrnrlat = 33.5
urcrnrlat = 37.2
llcrnrlon = -103.2
urcrnrlon = -94.0

# Initialize map
print 'Initializing map'
fig_map, ax_map = plt.subplots(1, figsize=(12, 6.75))
m = Basemap(projection='merc', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, 
    llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon, resolution='l', ax=ax_map)

# Draw patches over not Oklahoma
print 'Reading shapefile'
m.readshapefile('/Users/briangreene/Nextcloud/thermo/documentation/States/states_21basic/states',
 'states', drawbounds = False)
patch = []
patch_ok = []
for info, shape in zip(m.states_info, m.states):
	if info['STATE_NAME'] != 'Oklahoma':
		patch.append( Polygon(np.array(shape), True) )
	else:
		patch_ok.append( Polygon(np.array(shape), True))

xy = patch_ok[0].get_xy()


print 'Creating meshgrid for map'
xnew = np.arange(llcrnrlon, urcrnrlon, 0.1)
ynew = np.arange(llcrnrlat, urcrnrlat, 0.1)
xx, yy = np.meshgrid(xnew, ynew)




# # Convert to map coordinates
# xplot, yplot = m(xxyy[:, 0], xxyy[:, 1])

# # convert back to shape(-1, 2)
# xxplt, yyplt = np.meshgrid(xplot, yplot)
# xxyyplt = np.dstack((xxplt, yyplt))
# xxyyplt = xxyyplt.reshape(-1, 2)
# def to_left(v1, v2):
#     return (v1[0]*v2[1] - v1[1]*v2[0]) > 0



# plt.plot(xxyy[:,0], xxyy[:,1])


# ax_map.add_collection(PatchCollection(patch, facecolor= 'k', edgecolor='k', 
# 	linewidths=1., zorder=2))

# Grab current Mesonet data
print 'Grabbing mesonet data'
URL = 'http://www.mesonet.org/data/public/mesonet/current/current.csv.txt'
fd = urllib2.urlopen(URL)
data_long = fd.read()
data = data_long.split('\n')
lat_, lon_, yr_, mo_, da_, hr_, mi_, T_, Td_, Tchill_, wind_, dir_, Temps, \
	Dewpoint, Chill, Wind, Dir, lats, lons = ([] for i in range(19))
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
	Td_.append(data_short[11])
	Tchill_.append(data_short[13])
	wind_.append(data_short[16])
	dir_.append(data_short[15])

for i in np.arange(len(T_)):
	if T_[i] != ' ':
		Temps.append(float(T_[i]))
		Dewpoint.append(float(Td_[i]))
		#Chill.append(float(Tchill_[i]))
		#Wind.append(float(wind_[i]))
		#Dir.append(dir_[i]) 
		lats.append(lat_[i])
		lons.append(lon_[i])   

# transpose to basemap coordinates
print 'transposing to basemap coords'
lon, lat, xplot, yplot = m.makegrid(len(xnew), len(ynew), returnxy=True)
x1, y1 = m(lons, lats)

# xxyy = np.dstack((xplot, yplot))
# xxyy = xxyy.reshape(-1, 2)
# xxplot, yyplot = np.meshgrid(xxyy[:,0], xxyy[:,1])

points = np.array(zip(x1, y1))

t0 = time.time()
Tnew = interpolate.griddata(points, Dewpoint, (xplot, yplot), method='cubic')
print "Interpolate time: %s seconds" % (time.time() - t0)
# Create mask
mpath = Path(xy)
xxyyplot = np.dstack((xplot, yplot))
xxyyplot_flat = xxyyplot.reshape(-1, 2)
tstart = time.time()
mask_flat = np.invert(mpath.contains_points(xxyyplot_flat))
mask = mask_flat.reshape(xplot.shape)
print "Mask time: %s seconds" % (time.time() - tstart)

Tplot = ma.masked_where(mask, Tnew)

cfax = m.pcolormesh(xplot, yplot, Tplot, cmap=cmocean.cm.haline, vmin=5., vmax=25.)
cbar = m.colorbar(cfax)

for i in np.arange(len(Temps)):
	plotstring = '%d' % Dewpoint[i]
	txt= plt.text(x1[i], y1[i], plotstring, fontweight='bold', 
		ha='center', va='center')
	txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

m.drawcounties()
m.drawstates()

title_str = 'Oklahoma Mesonet 2m Td - %d%d%d - %d:%d CST' % (yr_[0], mo_[0], 
	da_[0], hr_[0], mi_[0])
ax_map.set_title(title_str)
cbar.ax.set_ylabel('Dewpoint Temperature ($^\circ$F)')

# ax_map.add_collection(PatchCollection(patch_ok, facecolor='b', edgecolor='r',
# 	linewidths=1., zorder=2))


plt.show()