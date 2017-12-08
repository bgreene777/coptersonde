import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import matplotlib.dates as mpdates
from matplotlib import cm
import cmocean
import numpy as np
import numpy.ma as ma
import csv
import math
import os
from glob import glob
from scipy.interpolate import interp1d
from datetime import datetime, timedelta
import Tkinter
from tkFileDialog import askdirectory
from metpy.constants import dry_air_gas_constant as Rd
from metpy.constants import water_heat_vaporization as Lv
from metpy.constants import dry_air_spec_heat_press as cp
import warnings

## Required pacakges & files: cmocean
## Change the following directories:
# ****** os filesep = '/' and windows filesp = '\\' ******
# ****** this version formatted for mac os ******

# Ignore warnings
warnings.filterwarnings("ignore", 
	".*converting a masked element to nan.*")
warnings.filterwarnings("ignore",
	".*invalid value encountered in greater_equal.*")
warnings.filterwarnings("ignore",
	".*invalid value encountered in true_divide.*")

# Constants
Rd = 1000. * Rd.magnitude
Lv = Lv.magnitude
cp = float(cp.magnitude)

# Select directory
root = Tkinter.Tk()
root.withdraw()
root.update()
initName = '/Users/briangreene/Nextcloud/thermo/data/PBL Transition/'
dataDirName = askdirectory(initialdir=initName)
root.destroy()

#dataDirName = '/Users/briangreene/Nextcloud/thermo/data/CLOUDMAP17/CSV/20170629/OSUF'
fnameArr = glob(os.path.join(dataDirName, "*.csv")) # string array

## Initialize
# Be smart enough to find number of heights automatically
nHtest = 1000
numfiles = len(fnameArr)
findHeights = np.full((nHtest, numfiles), np.nan)
count = 0
for f in fnameArr:
	H = []
	a = open(f)
	reader = csv.DictReader(a)
	for line in reader:
		H.append(float(line['AltAGL(m)']))

	a.close()
	size = len(H)
	findHeights[:size, count] = H
	count += 1

# Now go back and read in data for realz
nHeights = int(np.nanmax(findHeights) / 10.)
time = []
#lat = np.full((nHeights, numfiles), np.nan)
#lon = np.full((nHeights, numfiles), np.nan)
alt = np.full((nHeights, numfiles), np.nan)
p = np.full((nHeights, numfiles), np.nan)
T = np.full((nHeights, numfiles), np.nan)
Td = np.full((nHeights, numfiles), np.nan)
RH = np.full((nHeights, numfiles), np.nan)
w = np.full((nHeights, numfiles), np.nan)
theta = np.full((nHeights, numfiles), np.nan)
speed = np.full((nHeights, numfiles), np.nan)
direction = np.full((nHeights, numfiles), np.nan)
numHeightsArr = []

## Loop through directory and import data
count = 0
for fname in fnameArr:

    if fname.endswith('.csv'):
    	print '>>Reading file %s' % fname.split('/')[-1]
    else:
        print '>>>Wtf are you doing'
        print fname
        break


    # Initialize instance variables
    #latitude = []
    #longitude = []
    altitude = []
    pressure = []
    temp = []
    dew = []
    rh = []
    mix = []
    pottemp = []
    spd = []
    wind = []
    # Open file
    f = open(fname)
    reader = csv.DictReader(f)
    # Assign
    for line in reader:
    	#latitude.append(float(line['Lat']))
    	#longitude.append(float(line['Lon']))
    	altitude.append(float(line['AltAGL(m)']))
    	pressure.append(float(line['p(hPa)']))
    	temp.append(float(line['T(C)']))
    	dew.append(float(line['Td(C)']))
    	rh.append(float(line['RH(percent)']))
    	# want w in kg kg-1
    	mix.append(float(line['w(gKg-1)'])/1000.)
    	pottemp.append(float(line['Theta(K)']))
    	spd.append(float(line['Speed(ms-1)']))
    	wind.append(float(line['Dir(deg)']))

    f.close()

    # Assign to arrays
    sz = len(temp)
    #lat[:sz, count] = latitude
    #lon[:sz, count] = longitude
    alt[:sz, count] = altitude
    p[:sz, count] = pressure
    T[:sz, count] = temp
    Td[:sz, count] = dew
    RH[:sz, count] = rh
    w[:sz, count] = mix
    theta[:sz, count] = pottemp
    speed[:sz, count] = spd
    direction[:sz, count] = wind

    # Mask values
    #lat = ma.masked_invalid(lat)
    #lon = ma.masked_invalid(lon)
    alt = ma.masked_invalid(alt)
    p = ma.masked_invalid(p)
    T = ma.masked_invalid(T)
    Td = ma.masked_invalid(Td)
    RH = ma.masked_invalid(RH)
    w = ma.masked_invalid(w)
    theta = ma.masked_invalid(theta)
    speed = ma.masked_invalid(speed)
    direction = ma.masked_invalid(direction)

    time.append(datetime.strptime(fname.split('/')[-1][0:15], '%Y%m%d_%H%M%S'))
    numHeightsArr.append(sz)

    count += 1


## Set up Dimensions
time_list = [i.strftime('%H:%M:%S') for i in time]
title_today = time[0].strftime('%d %b %Y')
maxHeight = max(numHeightsArr)

numInterp = 1.
delta_t = 60. * numInterp
heights = 10. * np.arange(0, nHeights)

dt = mpdates.date2num(time)
dt_interp = mpdates.drange(time[0], time[-1], timedelta(seconds=delta_t))

## Calculate sensible and latent heat fluxes
dh = 10. # m

# time intervals
delta_t = np.array([(time[i] - time[i-1]).total_seconds() 
	for i in np.arange(1, len(dt))])

# initialize
delta_theta = np.full((nHeights, numfiles-1), np.nan)
delta_q = np.full((nHeights, numfiles-1), np.nan)
H0 = np.full((nHeights, numfiles), np.nan)
F0 = np.full((nHeights, numfiles), np.nan)
rho = np.full((nHeights, numfiles), np.nan)
q = np.full((nHeights, numfiles), np.nan)
H = np.full((nHeights, numfiles), np.nan)
F = np.full((nHeights, numfiles), np.nan)

# density and specific humidity
for j in np.arange(numfiles):
	for i in np.arange(nHeights):
		rho[i, j] = 100. * p[i, j] / (Rd * (T[i, j] + 273.15))
		q[i, j] = w[i, j] / (1. + w[i, j])

# change in theta and q in subsequent times at each level
for j in np.arange(1, numfiles):
	for i in np.arange(nHeights):
		delta_theta[i-1, j-1] = theta[i, j] - theta[i, j-1]
		delta_q[i-1, j-1] = q[i, j] - q[i, j-1]

# H and F at each level and time
for j in np.arange(1, numfiles):
	for i in np.arange(1, nHeights):
		H0[i, j] = cp * rho[i, j] * dh * delta_theta[i-1, j-1] / delta_t[j-1]
		F0[i, j] = Lv * rho[i, j] * dh * delta_q[i-1, j-1] / delta_t[j-1]

# sum terms from top down
for j in np.arange(0, numfiles):
	for i in np.arange(0, nHeights):
		H[i, j] = np.nansum(H0[i:, j])
		F[i, j] = np.nansum(F0[i:, j])

## Calculate static stability
dtdz = np.full((nHeights, numfiles), np.nan)
alt_windmax = np.full((1, numfiles), np.nan)
for j in np.arange(0, numfiles):
	alt_windmax[0, j] = alt[np.argmax(speed[:, j]), j]
	for i in np.arange(1, nHeights):
		dtdz[i, j] = (theta[i, j] - theta[i-1, j]) / dh


## Set up plotting parameters
# Alpha levels
alpha = np.linspace(0.1, 1, num=numfiles)

# Meshgrid for barbs
t_barbs = int(10. / numInterp)
xx, yy = np.meshgrid(dt_interp[0::t_barbs], heights[0::4])
xreal, yreal = np.meshgrid(dt, heights[0::4])

# Interpolate
def interpTime(t, tnew, z, data):
	data_interp = np.full((len(z),len(tnew)), np.nan)
	for i in np.arange(len(z)):
		f = interp1d(t,data[i,:])
		fnew = f(tnew)
		data_interp[i,:] = fnew

	data_interp = ma.masked_invalid(data_interp)
	return data_interp

T_interp = interpTime(dt, dt_interp, heights, T)
Td_interp = interpTime(dt, dt_interp, heights, Td)
RH_interp = interpTime(dt, dt_interp, heights, RH)
w_interp = interpTime(dt, dt_interp, heights, w)
theta_interp = interpTime(dt, dt_interp, heights, theta)
# need to interpolate components, then convert back if necessary
u = np.multiply(-speed, np.sin(direction * np.pi / 180.))
u_interp_kts = interpTime(dt, dt_interp, heights, u) * 1.94
v = np.multiply(-speed, np.cos(direction * np.pi / 180.))
v_interp_kts = interpTime(dt, dt_interp, heights, v) * 1.94
dtdz_interp = interpTime(dt, dt_interp, heights, dtdz)

f = interp1d(dt, alt_windmax)
alt_windmax_interp = f(dt_interp).reshape((-1,))

# Contour levels
Tlevels = np.arange(np.floor(np.nanmin(T_interp)),
	np.ceil(np.nanmax(T_interp)) + 1.0, 1.0)
Tdlevels = np.arange(np.floor(np.nanmin(Td_interp)),
	np.ceil(np.nanmax(Td_interp)) + 1.0, 1.0)
RHlevels = np.arange(np.floor(np.nanmin(RH_interp)),
	np.ceil(np.nanmax(RH_interp)) + 5.0, 5.0)
wlevels = np.arange(np.floor(np.nanmin(w_interp)),
	np.ceil(np.nanmax(w_interp)) + 0.5, 0.5)
thetalevels = np.arange(np.floor(np.nanmin(theta_interp)),
	np.ceil(np.nanmax(theta_interp)) + 1.0, 1.0)

# Check for existence of data to plot
if Tlevels.size == 0:
	isPlotT = 0
else:
	isPlotT = 1
if Tdlevels.size == 0:
	isPlotTd = 0
else:
	isPlotTd = 1
if RHlevels.size == 0:
	isPlotRH = 0
else:
	isPlotRH = 1
if wlevels.size == 0:
	isPlotw = 0
else:
	isPlotw = 1
if thetalevels.size == 0:
	isPlottheta = 0
else:
	isPlottheta = 1

## Plot
# Temperature
if isPlotT:
	fig1, ax1 = plt.subplots(1, figsize = (16,9))
	cfax1 = plt.pcolormesh(dt_interp, heights, T_interp,
		cmap=cmocean.cm.thermal)
	cax1 = plt.contour(dt_interp, heights, T_interp, levels=Tlevels, colors='k')
	plt.clabel(cax1, fontsize=9, inline=1, fmt='%3.1f')
	ax1.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
	ax1.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
	plt.xlabel('Time UTC')
	plt.ylabel('Altitude AGL (m)')
	plt.title('Temperature ' + title_today)
	cbar1 = fig1.colorbar(cfax1)
	cbar1.ax.set_ylabel('Temperature ($^\circ$C)')
	[plt.axvline(t, linestyle='--', color='k') for t in dt]

# Dewpoint
if isPlotTd:
	fig2, ax2 = plt.subplots(1, figsize=(16,9))
	cfax2 = plt.pcolormesh(dt_interp, heights, Td_interp, 
		cmap=cmocean.cm.haline)
	cax2 = plt.contour(dt_interp, heights, Td_interp, levels=Tdlevels, colors='k')
	plt.clabel(cax2, fontsize=9, inline=1, fmt='%3.1f')
	ax2.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
	ax2.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
	plt.xlabel('Time UTC')
	plt.ylabel('Altitude AGL (m)')
	plt.title('Dewpoint Temperature ' + title_today)
	cbar2 = fig2.colorbar(cfax2)
	cbar2.ax.set_ylabel('Dewpoint Temperature ($^\circ$C)')
	[plt.axvline(t, linestyle='--', color='k') for t in dt]

# RH
if isPlotRH:
	fig3, ax3 = plt.subplots(1, figsize=(16,9))
	cfax3 = plt.pcolormesh(dt_interp, heights, RH_interp,
		cmap=cmocean.cm.haline)
	cax3 = plt.contour(dt_interp, heights, RH_interp, levels=RHlevels, colors='k')
	plt.clabel(cax3, fontsize=9, inline=1, fmt='%g')
	ax3.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
	ax3.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
	plt.xlabel('Time UTC')
	plt.ylabel('Altitude AGL (m)')
	plt.title('Relative Humidity ' + title_today)
	cbar3 = fig3.colorbar(cfax3)
	cbar3.ax.set_ylabel('Relative Humidity (%)')
	[plt.axvline(t, linestyle='--', color='k') for t in dt]

# Potential Temperature
if isPlottheta:
	fig4, ax4 = plt.subplots(1, figsize=(16,9))
	cfax4 = plt.pcolormesh(dt_interp, heights, theta_interp, 
		cmap=cmocean.cm.thermal)
	cax4 = plt.contour(dt_interp, heights, theta_interp, levels=thetalevels, colors='k')
	plt.clabel(cax4, fontsize=9, inline=1, fmt='%g')
	ax4.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
	ax4.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
	plt.xlabel('Time UTC')
	plt.ylabel('Altitude AGL (m)')
	plt.title('Potential Temperature ' + title_today)
	cbar4 = fig4.colorbar(cfax4)
	cbar4.ax.set_ylabel('Potential Temperature (Kelvin)')
	[plt.axvline(t, linestyle='--', color='k') for t in dt]

# Mixing Ratio
if isPlotw:
	fig5, ax5 = plt.subplots(1, figsize=(16,9))
	cfax5 = plt.pcolormesh(dt_interp, heights, w_interp, 
		cmap=cmocean.cm.haline)
	cax5 = plt.contour(dt_interp, heights, w_interp, levels=wlevels, colors='k')
	plt.clabel(cax5, fontsize=9, inline=1, fmt='%3.1f')
	ax5.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
	ax5.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
	plt.xlabel('Time UTC')
	plt.ylabel('Altitude AGL (m)')
	plt.title('Mixing Ratio ' + title_today)
	cbar5 = fig5.colorbar(cfax5)
	cbar5.ax.set_ylabel('Mixing Ratio (g Kg-1)')
	[plt.axvline(t, linestyle='--', color='k') for t in dt]

# Wind Barbs
fig6, ax6 = plt.subplots(1, figsize=(16,9))
ax6.barbs(xx, yy, u_interp_kts[0::4,0::t_barbs], v_interp_kts[0::4,0::t_barbs])
#ax6.barbs(xreal, yreal, u[0::4,:], v[0::4,:])
ax6.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
ax6.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
plt.xlabel('Time UTC')
plt.ylabel('Altitude AGL (m)')
plt.title('Wind Speed (kts) and Direction ' + title_today)
[plt.axvline(t, linestyle='--', color='r') for t in dt]

# Sensible heat flux
fig7, ax7 = plt.subplots(1, figsize=(8,8))
plt.xlabel('Sensible Heat Flux (W m$^{-2}$)')
plt.ylabel('Altitude AGL (m)')
plt.title('Sensible Heat Flux Evolution')
[plt.plot(H[:, iplot], alt[:, iplot], label=time_list[iplot],
	color=(1, 0, 0, alpha[iplot])) 
	for iplot in np.arange(1, numfiles)]
ax7.legend()

# Latent heat flux
fig8, ax8 = plt.subplots(1, figsize=(8,8))
plt.xlabel('Latent Heat Flux (W m$^{-2}$)')
plt.ylabel('Altitude AGL (m)')
plt.title('Latent Heat Flux Evolution')
[plt.plot(F[:, iplot], alt[:, iplot], label=time_list[iplot],
	color=(0, 0, 1, alpha[iplot])) 
	for iplot in np.arange(1, numfiles)]
ax8.legend()

# Static Stability, Jet max altitude
fig9, ax9 = plt.subplots(1, figsize=(16,9))
cfax9 = plt.pcolormesh(dt_interp, heights, dtdz_interp, 
	cmap=cmocean.cm.balance, vmin=-0.25, vmax=0.25)
ax9.plot(dt_interp, alt_windmax_interp, linewidth=3, color='k')
ax9.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
ax9.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
plt.xlabel('Time UTC')
plt.ylabel('Altitude AGL (m)')
plt.title('dtheta/dz and wind max altitude ' + title_today)
cbar9 = fig9.colorbar(cfax9)
cbar9.ax.set_ylabel('dtheta/dz (K m$^{-1}$)')
[plt.axvline(t, linestyle='--', color='k') for t in dt]

plt.show(block=False)

## Save 
s = raw_input('>>Save images? y/n ')
while s != 'y' and s != 'n':
	s = raw_input('>>Save images? y/n ')

if s == 'y':
	saveNameBase = time[0].strftime('%Y%m%d') + '-' + dataDirName[-4:]
	dirName = '/Users/briangreene/Documents/Coptersonde/KAEFS/Figures/' + \
		time[0].strftime('%Y%m%d') + '/'

	#savePNG1 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'Temperature')
	savePNG1 = '%s%s_%s' % (dirName, saveNameBase, 'Temperature')
	fig1.savefig(savePNG1)
	print '>>Saving file %s' % savePNG1.split('/')[-1]
	#savePNG2 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'Dewpoint')
	savePNG2 = '%s%s_%s' % (dirName, saveNameBase, 'Dewpoint')
	fig2.savefig(savePNG2)
	print '>>Saving file %s' % savePNG2.split('/')[-1]
	#savePNG3 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'RH')
	savePNG3 = '%s%s_%s' % (dirName, saveNameBase, 'RH')
	fig3.savefig(savePNG3)
	print '>>Saving file %s' % savePNG3.split('/')[-1]
	#savePNG4 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'Theta')
	savePNG4 = '%s%s_%s' % (dirName, saveNameBase, 'Theta')
	fig4.savefig(savePNG4)
	print '>>Saving file %s' % savePNG4.split('/')[-1]
	#savePNG5 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'Mixing')
	savePNG5 = '%s%s_%s' % (dirName, saveNameBase, 'Mixing')
	fig5.savefig(savePNG5)
	print '>>Saving file %s' % savePNG5.split('/')[-1]
	#savePNG6 = '/Users/briangreene/Desktop/%s_%s' % (saveNameBase, 'Wind')
	savePNG6 = '%s%s_%s' % (dirName, saveNameBase, 'Wind')
	fig6.savefig(savePNG6)
	print '>>Saving file %s' % savePNG6.split('/')[-1]

	print '>>Figures Saved!'
elif s == 'n':
	print '>>Figures not saved. '
else:
	print '>>How did you get here??'

## Quit when ready
q = raw_input('>>Press enter to quit. ')
while q != '':
    q = raw_input('>>>Press enter to quit. ')

plt.close('all')
print 'Analysis Complete.'