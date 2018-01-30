from matplotlib import pyplot as plt
import matplotlib.dates as mpdates
from matplotlib import cm
import cmocean
import numpy as np
import numpy.ma as ma
import netCDF4
import math
import os
import sys
from glob import glob
from scipy.interpolate import interp1d, UnivariateSpline
from datetime import datetime, timedelta
from metpy.constants import dry_air_gas_constant as Rd
from metpy.constants import water_heat_vaporization as Lv
from metpy.constants import dry_air_spec_heat_press as cp
import warnings
import wxtools

'''
Reads nc file of combined profiles
Prompts user for date (YYYYMMDD) and location of profiles
Brian Greene
University of Oklahoma
Last edit: 29 January 2018
'''

## Required pacakges & files: cmocean

# Ignore warnings
warnings.filterwarnings("ignore", 
	".*converting a masked element to nan.*")
np.seterr(invalid='ignore')
warnings.filterwarnings("ignore",
	".*invalid value encountered in greater_equal.*")
warnings.filterwarnings("ignore",
	".*invalid value encountered in true_divide.*")

# Constants
Rd = 1000. * Rd.magnitude
Lv = Lv.magnitude
cp = float(cp.magnitude)

# Directory from input argument
#inputs = sys.argv
procDate = raw_input('>>Enter date of profiles (YYYYMMDD): ')
procLoc = raw_input('>>Enter location: ')

# First check if composite nc file exists
user, isMac = wxtools.getUser()
dataDirName = os.sep + os.path.join('Users', user, 'Nextcloud',
	'thermo', 'data', 'PBL Transition', procDate, procLoc, 'pp_nc')

# Filename based off user input
fname = glob(os.path.join(dataDirName, "*all.nc")) # string array

if len(fname) == 0:
	# if this file does not exist, run csvdir_to_ncdir and collect_nc
	print '>>nc files do not exist. Initializing csvdir_to_ncdir '
	dataDircsv = os.path.join(dataDirName.rsplit(os.sep, 1)[0], 'pp_csv')
	wxtools.csvdir_to_ncdir(dataDircsv + os.sep)
	print '>>Initializing collect_nc '
	wxtools.collect_nc(procDate, procLoc)
	fname = glob(os.path.join(dataDirName, "*all.nc"))[0]
else:
	fname = fname[0]

# Import nc file
print '>>Reading file: %s' % fname.split(os.sep)[-1]
df = netCDF4.Dataset(fname, 'r')
# Optional data dictionary for irregular datasets
# dataDic = {}
# varNames = df.variables.keys()
# for key in varNames:
# 	dataDic.update({key:[]})
# 	dataDic[key].append(df.variables[key][:, :])

# Assign datasets
nHeights = len(df.variables['level'])
maxHeight = df.max_alt_agl_m
t = df.variables['time'][:]
numfiles = len(t)
time = []
[time.append(datetime.fromtimestamp(t[i])) for i in range(numfiles)]
alt = df.variables['alt_agl'][:, :].filled(np.nan)
p = df.variables['pressure'][:, :].filled(np.nan)
T = df.variables['Temperature'][:, :].filled(np.nan)
Td = df.variables['Dewpoint'][:, :].filled(np.nan)
RH = df.variables['RH'][:, :].filled(np.nan)
w = df.variables['Mixing'][:, :].filled(np.nan)
theta = df.variables['Theta'][:, :].filled(np.nan)
speed = df.variables['Speed'][:, :].filled(np.nan)
direction = df.variables['Dir'][:, :].filled(np.nan)

## Set up Dimensions
time_list = [i.strftime('%H:%M:%S') for i in time]
title_today = time[0].strftime('%d %b %Y')
#maxHeight = max(numHeightsArr)

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
	alt_windmax[0, j] = alt[np.nanargmax(speed[:, j]), j]
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
	return ma.masked_invalid(data_interp)

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

# # Sensible heat flux
# fig7, ax7 = plt.subplots(1, figsize=(8,8))
# plt.xlabel('Sensible Heat Flux (W m$^{-2}$)')
# plt.ylabel('Altitude AGL (m)')
# plt.title('Sensible Heat Flux Evolution')
# [plt.plot(H[:, iplot], alt[:, iplot], label=time_list[iplot],
# 	color=(1, 0, 0, alpha[iplot])) 
# 	for iplot in np.arange(1, numfiles)]
# ax7.legend()

# # Latent heat flux
# fig8, ax8 = plt.subplots(1, figsize=(8,8))
# plt.xlabel('Latent Heat Flux (W m$^{-2}$)')
# plt.ylabel('Altitude AGL (m)')
# plt.title('Latent Heat Flux Evolution')
# [plt.plot(F[:, iplot], alt[:, iplot], label=time_list[iplot],
# 	color=(0, 0, 1, alpha[iplot])) 
# 	for iplot in np.arange(1, numfiles)]
# ax8.legend()

# Static Stability, Jet max altitude
vmax = np.round(np.max(np.abs(dtdz_interp)), 2)
fig9, ax9 = plt.subplots(1, figsize=(16,9))
cfax9 = plt.pcolormesh(dt_interp, heights, dtdz_interp, 
	cmap=cmocean.cm.balance, vmin=-vmax, vmax=vmax)
ax9.plot(dt_interp, alt_windmax_interp, linewidth=3, color='k')
ax9.xaxis.set_major_locator(mpdates.MinuteLocator(interval=30))
ax9.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
plt.xlabel('Time UTC')
plt.ylabel('Altitude AGL (m)')
plt.title('dtheta/dz and wind max altitude ' + title_today)
cbar9 = fig9.colorbar(cfax9)
cbar9.ax.set_ylabel('dtheta/dz (K m$^{-1}$)')
[plt.axvline(t, linestyle='--', color='k') for t in dt]

#plt.show()
# plt.draw()
# plt.pause(0.01)

# Save 
s = raw_input('>>Save images? y/n ')
while s != 'y' and s != 'n':
	s = raw_input('>>Save images? y/n ')

if s == 'y':
	# Save name
	saveNameBase = df.date_str + '_' + df.location_short
	saveDirName = os.sep + os.path.join('Users', os.getlogin(), 'Desktop', 
		saveNameBase)
	# Check for existence of folder
	if not os.path.exists(saveDirName):
		os.mkdir(saveDirName)

	savePNG1 = os.path.join(saveDirName, saveNameBase + '_Temperature.png')
	fig1.savefig(savePNG1)
	print '>>Saving file %s' % savePNG1.split(os.sep)[-1]
	savePNG2 = os.path.join(saveDirName, saveNameBase + '_Dewpoint.png')
	fig2.savefig(savePNG2)
	print '>>Saving file %s' % savePNG2.split(os.sep)[-1]
	savePNG3 = os.path.join(saveDirName, saveNameBase + '_RH.png')
	fig3.savefig(savePNG3)
	print '>>Saving file %s' % savePNG3.split(os.sep)[-1]
	savePNG4 = os.path.join(saveDirName, saveNameBase + '_Theta.png')
	fig4.savefig(savePNG4)
	print '>>Saving file %s' % savePNG4.split(os.sep)[-1]
	savePNG5 = os.path.join(saveDirName, saveNameBase + '_Mixing.png')
	fig5.savefig(savePNG5)
	print '>>Saving file %s' % savePNG5.split(os.sep)[-1]
	savePNG6 = os.path.join(saveDirName, saveNameBase + '_Wind.png')
	fig6.savefig(savePNG6)
	print '>>Saving file %s' % savePNG6.split(os.sep)[-1]
	savePNG9 = os.path.join(saveDirName, saveNameBase + '_Static_Stability.png')
	fig9.savefig(savePNG9)
	print '>>Saving file %s' % savePNG9.split(os.sep)[-1]

	print '>>Figures Saved!'
elif s == 'n':
	print '>>Figures not saved. '
else:
	print '>>How did you get here??'

sh = raw_input('>>Show images? y/n ')
while sh != 'y' and sh != 'n':
	sh = raw_input('>>Show images? y/n ')

if sh == 'y':
	print 'Close figures to finish. '
	plt.show()

## Quit when ready
#q = raw_input('>>Press enter to quit. ')
# q = []
# while q != '':
# 	plt.pause(2)
# 	q = raw_input('>>>Press enter to quit. ')

# plt.close('all')
print 'Analysis Complete.'