import numpy as np
import numpy.ma as ma
import os
from glob import glob
import csv
import matplotlib.pyplot as plt
import matplotlib.dates as mpdates
import matplotlib.gridspec as gridspec
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
from metpy.plots import Hodograph, SkewT
import metpy.calc as mcalc
from metpy.units import units
from metpy.constants import dry_air_gas_constant as Rd
import geopy.distance
import cmocean
import math
from datetime import datetime, timedelta
import time
import urllib2
import pytz
import netCDF4

#############
## Version ##
#############
'''
Updated 21 December 2017
Brian Greene
University of Oklahoma
Various subroutines for use with OU Coptersonde to atuomate calculations
'''
#################
## Directories ##
#################

myNextcloud = '/Users/briangreene/Nextcloud/thermo/'
# location of raw .csv and .pos files
dataDirName = myNextcloud + 'data/'
# location of mesonet location csv
mesocsv = myNextcloud + 'documentation/Mesonet/geoinfo.csv'
# location of Logos.png
logoName = myNextcloud + 'documentation/Logos.png'

######################################
## Thermodynamics and Cloud Physics ##
######################################

def e_sl(T_C):
	'''
	Input: Temperature (Celsius)
	Goff-Gratch solution to Clausius-Clapeyron for liquid water
	Returns es (hPa)
	'''
	T = T_C + 273.15
	log_es = -7.90298 * ((373.15 / T) - 1.) + 5.02808 * np.log10(373.15/T)\
		- 1.3816 * (10 ** -7) * (10**(11.344 * (1. - T / 373.15)) - 1.)\
		+ 8.1328 * (10 ** -3) * (10**(-3.49149 * ((373.15/T) - 1.)) - 1.)\
		+ np.log10(1013.25)
	return 10**log_es

def e_si(T_C):
	'''
	Input: T (Celsius)
	Goff-Gratch solution to Clausius-Clapeyron for ice
	Returns ei (hPa)
	'''
	T = T_C + 273.15
	log_ei = -9.09718 * ((273.15 / T) - 1.) - 3.56654 * np.log10(273.15 / T)\
		+ 0.876793 * (1. - (T / 273.15)) + np.log10(6.11)
	return 10**log_ei

###################
## UAV Functions ##
###################

def csvread_raw(coptersondefilename):
	'''
	Input: filepath of coptersonde csv
	Returns numpy array of select coptersonde csv data

	Headers = Time, Lat, Long, Altitude, Abs Pressure, Roll, Pitch, Yaw, 
	Humidity_1, Humidity_2, Humidity_3, Humidity_4, Temp_1, Temp_2, Temp_3,
	Temp_4, Gnd Speed, Vert Speed, AccelX, AccelY
	'''
	f = open(coptersondefilename)
	reader = csv.DictReader(f)

	# Initialize
	time_, lat_, lon_, alt_, p_, roll_, pitch_, yaw_, RH1_, RH2_, RH3_, RH4_,\
		RHT1_, RHT2_, RHT3_, RHT4_, T1_, T2_, T3_, T4_, gspd_, vspd_, ax_, ay_ \
		= ([] for i in range(24))

	# Assign
	for line in reader:
	    time_.append(float(line['Time']))
	    lat_.append(float(line['Lat']))
	    lon_.append(float(line['Long']))
	    alt_.append(float(line['Altitude']))
	    p_.append(float(line['Abs Pressure']))
	    roll_.append(float(line['Roll']))
	    pitch_.append(float(line['Pitch']))
	    yaw_.append(float(line['Yaw']))
	    RH1_.append(float(line['Humidity_1']))
	    RH2_.append(float(line['Humidity_2']))
	    RH3_.append(float(line['Humidity_3']))
	    RH4_.append(float(line['Humidity_4']))
	    RHT1_.append(float(line['RH_1 Temp']))
	    RHT2_.append(float(line['RH_2 Temp']))
	    RHT3_.append(float(line['RH_3 Temp']))
	    RHT4_.append(float(line['RH_4 Temp']))
	    T1_.append(float(line['Temp_1']))
	    T2_.append(float(line['Temp_2']))
	    T3_.append(float(line['Temp_3']))
	    T4_.append(float(line['Temp_4']))
	    gspd_.append(float(line['Gnd Speed']))
	    vspd_.append(float(line['Vert Speed']))
	    ax_.append(float(line['AccelX']))
	    ay_.append(float(line['AccelY']))

	f.close()

	return np.array([time_, lat_, lon_, alt_, p_, roll_, pitch_, yaw_, 
		RH1_, RH2_, RH3_, RH4_, RHT1_, RHT2_, RHT3_, RHT4_, T1_, T2_, T3_, T4_, 
		gspd_, vspd_, ax_, ay_])

def dgpsread(ppkfilename):
	'''
	Input: filepath to .pos ppk file
	Output: numpy array of date_dgps, time_dgps, lat_dgps, long_dgps, alt_dgps, 
	alt_dgps_agl

	Be sure to set date_dgps and time_dgps to list
	Convert rest with .astype(np.float)
	'''
	fd = open(ppkfilename, 'r')
	for i in range(25):
		fd.next()

	reader = csv.DictReader(fd, delimiter=' ')
	reader.fieldnames = ['Date', 'Time', '', '', 'Lat', '', 'Lon', '', 'Alt2', 
    'Alt', '', '', 'Q', '', '', 'ns', '', '', 'sdn', '', '', 'sde', '', '',
    'sdu', '', '', 'sdne', '', '', 'sdeu', '', '', 'sdun', '', 'age', '', '',
    '', 'ratio']

	date_dgps_ = []
	time_dgps_ = []
	lat_dgps_ = []
	lon_dgps_ = []
	alt_dgps_ = []

	for line in reader:
	    date_dgps_.append(line['Date'])
	    time_dgps_.append(line['Time'])
	    lat_dgps_.append(float(line['Lat']))
	    lon_dgps_.append(float(line['Lon']))
	    # Altitude has different spaces once reaching 1000 ft
	    if not line['Alt']:
	        alt_dgps_.append(float(line['Alt2']))
	    else:
	        alt_dgps_.append(float(line['Alt']))

	fd.close()

	alt_dgps_agl_ = [a - alt_dgps_[0] for a in alt_dgps_]

	return np.array([date_dgps_, time_dgps_, lat_dgps_, lon_dgps_, alt_dgps_,
		alt_dgps_agl_])

def csvread_copter(coptercsv):
	'''
	Input: filepath of post-processed coptersonde csv
	Returns numpy array of coptersonde csv data using numpy loadtxt

	Headers = Lat, Lon, AltAGL(m), p(hPa), T(C), Td(C), RH(percent), w(gKg-1), 
	Theta(K), Speed(ms-1), Dir(deg)
	'''
	return np.loadtxt(coptercsv, delimiter=',', skiprows=1)

def csv_to_nc(filepath):
	'''
	Converts individual csv RAOB output to netcdf format
	Input filepath for csv
	Outputs new file in same directory
	'''
	# Load data
	data = csvread_copter(filepath)
	lat_ = data[0, 0]
	lon_ = data[0, 1]
	alt_ = data[:, 2]
	p_ = data[:, 3]
	T_ = data[:, 4]
	Td_ = data[:, 5]
	RH_ = data[:, 6]
	w_ = data[:, 7]
	Theta_ = data[:, 8]
	Speed_ = data[:, 9]
	Dir_ = data[:, 10]

	# Number of levels
	nHeights = len(alt_)
	# Grab filename, change .csv to .nc
	fname = filepath.split('/')[-1].split('.')[0] + '.nc'
	# Grab date and time from filename
	date_str, time_str = fname.split('-')[0].split('_')
	dt_ = datetime.strptime(fname.split('-')[0], '%Y%m%d_%H%M%S')
	timestamp_ = time.mktime(dt_.timetuple())
	# Find nearest mesonet site
	meso_ = findSite(lat_, lon_)
	elev_ = float(meso_.split('/')[2])
	# Covert AGL to MSL
	altMSL_ = alt_ + elev_

	# Initialize nc file
	rootgrp = netCDF4.Dataset(fname, 'w', format='NETCDF4')
	# Initialize dimensions - level, lat, lon
	level = rootgrp.createDimension('level', None)
	lat = rootgrp.createDimension('lat', 1)
	lon = rootgrp.createDimension('lon', 1)
	# Initialize variables
	# Dimensional info
	latitude = rootgrp.createVariable('lat', 'f4', ('lat',))
	longitude = rootgrp.createVariable('lon', 'f4', ('lon',))
	levels = rootgrp.createVariable('level', 'i4', ('level',))
	# Data
	alt_agl = rootgrp.createVariable('alt_agl', 'f4', ('level',))
	alt_msl = rootgrp.createVariable('alt_msl', 'f4', ('level',))
	p = rootgrp.createVariable('pressure', 'f4', ('level',))
	T = rootgrp.createVariable('Temperature', 'f4', ('level',))
	Td = rootgrp.createVariable('Dewpoint', 'f4', ('level',))
	RH = rootgrp.createVariable('RH', 'f4', ('level',))
	w = rootgrp.createVariable('Mixing', 'f4', ('level',))
	Theta = rootgrp.createVariable('Theta', 'f4', ('level',))
	Speed = rootgrp.createVariable('Speed', 'f4', ('level',))
	Dir = rootgrp.createVariable('Dir', 'f4', ('level',))
	# Global attributes
	rootgrp.decription = 'OU Coptersonde vertical profile'
	rootgrp.location_short = meso_.split('/')[0]
	rootgrp.location_long = meso_.split('/')[1]
	rootgrp.elevation_MSL_m = elev_
	rootgrp.date_str = date_str
	rootgrp.time_str = time_str
	rootgrp.timestamp = timestamp_
	rootgrp.max_alt_agl_m = np.nanmax(alt_)
	# Variable attributes
	alt_agl.name_long = 'Altitude above ground level'
	alt_agl.units = 'm'
	alt_msl.name_long = 'Altitude above sea level'
	alt_msl.units = 'm'
	p.name_long = 'Pressure'
	p.units = 'hPa'
	T.name_long = 'Temperature'
	T.units = 'C'
	Td.name_long = 'Dewpoint Temperature'
	Td.units = 'C'
	RH.name_long = 'Relative Humidity'
	RH.units = 'percent'
	w.name_long = 'Water Vapor Mixing Ratio'
	w.units = 'g kg-1'
	Theta.name_long = 'Potential Temperature'
	Theta.units = 'K'
	Speed.name_long = 'Wind Speed'
	Speed.units = 'm s-1'
	Dir.name_long = 'Wind Direction'
	Dir.units = 'degrees'
	# Assign values
	latitude[0] = lat_
	longitude[0] = lon_
	alt_agl[:nHeights] = alt_
	alt_msl[:nHeights] = altMSL_
	p[:nHeights] = p_
	T[:nHeights] = T_
	Td[:nHeights] = Td_
	RH[:nHeights] = RH_
	w[:nHeights] = w_
	Theta[:nHeights] = Theta_
	Speed[:nHeights] = Speed_
	Dir[:nHeights] = Dir_
	# Assign to unlimited dimension variable
	levels = np.arange(nHeights)

	# Close File
	rootgrp.close()

	print '>>%s created successfully!' % fname
	return


def findLatestDir(dirname):
	'''
	Input data directory for project
	Returns filepath of most recently created directory
	'''
	fnameArr = glob(os.path.join(dirname, '*/*/'))
	dirCreatedArr = []
	[dirCreatedArr.append(float(d.split('/')[-3])) for d in fnameArr]
	return fnameArr[max(enumerate(dirCreatedArr))[0]] + 'Raw/'

def findLatestCSV(dirname):
	'''
	Input data directory where raw csv files are saved
	Returns filepath of most recently created file
	'''
	latestDir = findLatestDir(dirname)
	fnameArr = glob(os.path.join(latestDir, '*.csv'))
	csvCreatedArr = []
	for s in fnameArr:
		csvCreatedArr.append(datetime.strptime(
			s.split('/')[-1][29:38], '%Hh%Mm%Ss'))
	return fnameArr[max(enumerate(csvCreatedArr))[0]]

def findDGPSfile(csvfilename):
	'''
	Input csv filepath, will read last datetime of csv and determine closest
	.pos file based on filename
	'''
	# read times from CSV
	dat = csvread_raw(csvfilename)
	timeCreateCSV = dat[0][-1]
	# Time created for all .pos files, based on file name
	dataDirName = csvfilename.rsplit('/', 1)[0]
	fnameArr = glob(os.path.join(dataDirName, "*.pos")) # string array
	posCreatedArr = []
	for s in fnameArr:
		idatetime = datetime.strptime(s.split('/')[-1][9:29],
			'%Y-%m-%d_%Hh%Mm%Ss')
		posCreatedArr.append(float(idatetime.strftime('%s')))
	posCreatedArr = np.array(posCreatedArr)
	# Find closest one

	return fnameArr[np.argmin(abs(posCreatedArr - timeCreateCSV))]

def findSite(cop_lat, cop_lon, fmeso=mesocsv):
	'''
	Inputs: copter latitude, copter longitude
	Locates nearest mesonet station using geopy vincenty distance
	Returns Site name, identifier, and elevation in form: SITE/longname/elev
	'''
	f = open(fmeso)
	reader = csv.DictReader(f)
	mesoLat, mesoLon, mesoName, mesoLong, mesoElev = ([] for i in range(5))
	for line in reader:
	    mesoLat.append(float(line['nlat']))
	    mesoLon.append(float(line['elon']))
	    mesoName.append(line['stid'])
	    mesoLong.append(line['name'])
	    mesoElev.append(line['elev'])
	mesoLat = np.array(mesoLat)
	mesoLon = np.array(mesoLon)
	mesoElev = np.array(mesoElev)
	# calculate distances to mesonet sites
	distances = np.array([geopy.distance.vincenty((cop_lat, cop_lon), (iPos)) \
	    for iPos in zip(mesoLat, mesoLon)])

	iMeso = distances.argmin()

	site = mesoName[iMeso]
	print 'Site identified: %s' % site
	return site + '/' + mesoLong[iMeso] + '/' + mesoElev[iMeso]

def findStart(vertSpdRaw, altRaw, timeRaw):
	'''
	Inputs: vertical speed, altitude, time
	Determines time to begin averageing. Finds time copter begins ascent after
	hovering and then selects time 5 points prior.
	Returns datetime format of beginning time
	'''
	# only look at first half of flight while hovering (below 15 m)
	i = np.argmax(altRaw)
	ii = np.where((altRaw[:i] > 8.) & (altRaw[:i] < 15.))
	althov = altRaw[ii]
	vspdhov = vertSpdRaw[ii]
	timehov = timeRaw[ii]
	# cut in half to only look at hovering -> ascent portion
	althov = althov[len(althov)/2:]
	vspdhov = vspdhov[len(vspdhov)/2:]
	timehov = timehov[len(timehov)/2:]
	iii = np.where(vspdhov[len(vspdhov)/2:] > 1.)[0][0] - 5
	return mpdates.num2date(timehov[iii])

def findEnd(vertSpdRaw, altRaw, timeRaw):
	'''
	Inputs: vertical speed, altitude, time
	Determines time to end averaging during descent. Finds time when copter's
	vertical speed changes from negative to positive.
	Returns datetime format of end time
	'''
	# only look at second half of flight
	i = np.argmax(altRaw)
	i += 100
	ii = np.squeeze(np.where(vertSpdRaw[i:] > 0.))[0]
	return mpdates.num2date(timeRaw[i:][ii])

def check_internet():
	try:
		urllib2.urlopen('http://129.15.109.106', timeout=1)
		return True
	except urllib2.URLError as err:
		return False

def getMesoData(procYear, procMonth, procDay, procStation):
	# first check for internet connection
	iswifi = check_internet()
	if iswifi:
		baseURL = 'http://mesonet.org/index.php/dataMdfMts/dataController/getFile/'
		URL = baseURL + '%4.4d%2.2d%2.2d%s/mts/TEXT/' % (procYear, procMonth, 
			procDay, procStation.lower())
		fd = urllib2.urlopen(URL)
		data_long = fd.read()
		data = data_long.split('\n')
		time_, RH_, T2m_, T9m_, speed_, direction_, p_, srad_ = ([] for i in range(8))

		for i in np.arange(3, len(data)-1):
			data_short = data[i].split()
			time_.append(float(data_short[2]))
			RH_.append(float(data_short[3]))
			T2m_.append(float(data_short[4]))
			T9m_.append(float(data_short[14]))
			speed_.append(float(data_short[5]))
			direction_.append(float(data_short[7]))
			p_.append(float(data_short[12]))
			srad_.append(float(data_short[13]))

		speed_kts = [i * 1.94 for i in speed_]
		u,v = mcalc.get_wind_components(speed_kts*units.kts, direction_ * units.deg)
		u = u.to(units.kts) / units.kts
		v = v.to(units.kts) / units.kts

		return np.array([time_, RH_, T2m_, T9m_, u, v, p_, srad_])

	else:
		return np.array([])

def findClosestMesoTime(timeCopter):
	'''
	timeCopter: datetime object of time takeoff or time land
	Returns index of nearest mesonet 5-minute observation
	'''
	mesoTimes = np.linspace(0, 1435, num=1440/5)
	baseDay = datetime(timeCopter.year, timeCopter.month, timeCopter.day)
	baseDay = baseDay.replace(tzinfo=pytz.UTC)
	timeCopter = timeCopter.replace(tzinfo=pytz.UTC)
	numMin = (timeCopter - baseDay).total_seconds() / 60.
	return np.argmin(abs(mesoTimes - numMin))

def uavCAPE(Tenv, Tprof, pressure):
	'''
	Input: Environmental Temperature, Profile Temperature, pressure
	Calculates surface-based cape based on following equation:
	CAPE = Rd * integral(T'(p) - T(p)) dlnp
	use "bottom Riemann sum":
	uavCAPE = Rd * sum(T'(p) - T(p)) * delta_p/p
	Returns SBCAPE for extent of profile
	'''
	dlnp = []
	for i in range(len(Tenv) - 1):
	    dlnp.append(np.abs(pressure[i+1] - pressure[i]) / pressure[i])

	dT = Tprof.to('kelvin') - Tenv.to('kelvin')
	dCAPE = dT[:-1] * dlnp
	return Rd.to('joule/kilogram/kelvin') * np.nansum(dCAPE) * dCAPE.units

def interpTime(t, tnew, z, data):
	'''
	Inputs: Time from data, interpolated time, selected altitudes, parameter to
	be displayed
	Interpolates data in time for each height for multiple UAV profiles
	Returns numpy array of interpolated values
	'''
	data_interp = np.full((len(z), len(tnew)), np.nan)
	for i in np.arange(len(z)):
		f = interp1d(t, data[i, :])
		fnew = f(tnew)
		data_interp[i, :] = fnew

	data_interp = ma.masked_invalid(data_interp)
	return data_interp

def parcelUAV(T, Td, p):
	'''
	Inputs: temperature, dewpoint, and pressure
	Returns: lcl pressure, lcl temperature, isbelowlcl flag, profile temp
	'''
	lclpres, lcltemp = mcalc.lcl(p[0] * units.mbar, 
        T[0] * units.degC, Td[0] * units.degC)
	print 'LCL Pressure: %5.2f %s' % (lclpres.magnitude, lclpres.units)
	print 'LCL Temperature: %5.2f %s' % (lcltemp.magnitude, lcltemp.units)

	# parcel profile
	# determine if there are points sampled above lcl
	ilcl = np.squeeze(np.where((p * units.mbar) <= lclpres))
	# if not, entire profile dry adiabatic
	if ilcl.size == 0:
	    prof = mcalc.dry_lapse(p * units.mbar, T[0] * units.degC).to('degC')
	    isbelowlcl = 1
	# if there are, need to concat dry & moist profile ascents
	else:
	    ilcl = ilcl[0]
	    prof_dry = mcalc.dry_lapse(p[:ilcl] * units.mbar,
	        T[0] * units.degC).to('degC')
	    prof_moist = mcalc.moist_lapse(p[ilcl:] * units.mbar,
	        prof_dry[-1]).to('degC')
	    prof = np.concatenate((prof_dry, prof_moist)) * units.degC
	    isbelowlcl = 0

	return lclpres, lcltemp, isbelowlcl, prof


def plotUAVskewT(filenamecsv, filename2=''):
	'''
	Input filepath of post-processed uav data
	Optionally enter filename of 2nd csv

	Outputs Skew-T log-p plot of UAV data, includes hodograph and some
	convective parameters
	'''
	# File 1
	copdata = csvread_copter(filenamecsv)
	lat = copdata[0]
	lon = copdata[1]
	alt = copdata[2]
	pressure = copdata[3]
	temperature = copdata[4]
	dewpoint = copdata[5]
	speed = copdata[9]
	speed_kts = speed * 1.94
	direction = copdata[10]
	site = findSite(lat[0], lon[0])
	sitename, sitelong = site.split('/')
	fname = filenamecsv.split('/')[-1]
	timeTakeoff = datetime.strptime(fname[:15], '%Y%m%d_%H%M%S')
	copterNum = fname[-10]

	u,v = mcalc.get_wind_components(speed_kts*units.kts, direction * units.deg)
	u = u.to(units.kts)
	v = v.to(units.kts)
	# Wind shear
	bulkshear = speed_kts[-3] - speed_kts[0]
	print '0-%d m Bulk Shear: %.0f kts' % (alt[-3], bulkshear)
	if np.isnan(dewpoint).all():
		moist = 0
	else: 
		moist = 1

	# File 2
	if len(filename2) > 0:
		copdata2 = csvread_copter(filename2)
		p2 = copdata2[3]
		T2 = copdata2[4]
		Td2 = copdata2[5]
		spd2 = copdata2[9]
		spd_kts2 = spd2 * 1.94
		dir2 = copdata2[10]
		u2, v2 = mcalc.get_wind_components(spd_kts2*units.kts, 
			dir2 * units.deg)
		u2 = u2.to(units.kts)
		v2 = v2.to(units.kts)
	else:
		p2, T2, Td2, Td2, spd2, spd_kts2, dir2, u2, v2 = ([] for i in range(9))



	print 'Plotting...'
	fignum = plt.figure(figsize=(12,9))
	gs = gridspec.GridSpec(4, 4)
	skew = SkewT(fignum, rotation=20, subplot=gs[:, :2])

	skew.plot(pressure, temperature, 'r', linewidth = 2)
	skew.plot(pressure, dewpoint, 'g', linewidth = 2)
	skew.plot_barbs(pressure[0::4], u[0::4], v[0::4], x_clip_radius = 0.12, \
	    y_clip_radius = 0.12)

	# Plot file 2
	if len(filename2) > 0:
		skew.plot(p2, T2, 'r', linewidth=2, alpha=0.5)
		skew.plot(p2, Td2, 'g', linewidth=2, alpha=0.5)
		skew.plot_barbs(p2[0::4], u2[0::4], v2[0::4], barbcolor='r')

	# Plot convective parameters
	if moist:
		plcl, Tlcl, isbelowlcl, profile = parcelUAV(temperature,
			dewpoint, pressure)
		SBCAPE = uavCAPE(temperature * units.degC, profile,
			pressure * units.hPa)
		skew.plot(plcl, Tlcl, 'ko', markerfacecolor='black')
		skew.plot(pressure, profile, 'k', linewidth=2)
	else:
		isbelowlcl = 0

	# set up plot limits and labels - use LCL as max if higher than profile
	if moist:
		xmin = math.floor(np.nanmin(dewpoint))
	else:
	    xmin = math.floor(np.nanmin(temperature))
	xmax = math.floor(np.nanmax(temperature)) + 20
	if isbelowlcl:
	    ymin = round((plcl / units.mbar), -1) - 10
	else:
	    ymin = round(np.nanmin(pressure),-1) - 10
	    
	ymax = round(np.nanmax(pressure),-1) + 10

	skew.ax.set_ylim(ymax, ymin)
	skew.ax.set_xlim(xmin, xmax)
	skew.ax.set_yticks(np.arange(ymin, ymax+10, 10))

	skew.ax.set_xlabel('Temperature ($^\circ$C)')
	skew.ax.set_ylabel('Pressure (hPa)')
	titleName = 'Coptersonde-%s %s UTC - %s' % (copterNum, 
	    timeTakeoff.strftime('%d-%b-%Y %H:%M:%S'), sitename)
	skew.ax.set_title(titleName)

	skew.plot_dry_adiabats(linewidth=0.75)
	skew.plot_moist_adiabats(linewidth=0.75)
	skew.plot_mixing_lines(linewidth=0.75)

	# Hodograph
	ax_hod = fignum.add_subplot(gs[:2,2:])
	#gs.tight_layout(fig5)
	if np.nanmax(speed_kts) > 18:
	    comprange = 35
	else:
	    comprange = 20

	h = Hodograph(ax_hod, component_range=comprange)
	h.add_grid(increment=5)
	h.plot_colormapped(u, v, pressure, cmap=cmocean.cm.deep_r)
	ax_hod.set_title('Hodograph (kts)')
	ax_hod.yaxis.set_ticklabels([])
	#ax_hod.set_xlabel('Wind Speed (kts)')

	# Map - Oklahoma
	llcrnrlat = 33.6
	urcrnrlat = 37.2
	llcrnrlon = -103.2
	urcrnrlon = -94.2
	ax_map = fignum.add_subplot(gs[2, 2:])

	m = Basemap(projection='merc', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, 
	    llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon, lat_ts=20, resolution='l',
	    ax=ax_map)

	print 'Basemap...'
	m.drawcounties()
	m.drawstates()
	x,y = m(lon[0], lat[0])
	plt.plot(x,y,'b.')
	plt.text(x+40000, y-5000, sitelong, 
		bbox=dict(facecolor='yellow', alpha=0.5))

	if moist:
	    # Convective parameter values
	    ax_data = fignum.add_subplot(gs[3, 2])
	    plt.axis('off')
	    datastr = 'LCL = %.0f hPa\nSBCAPE = %.0f J kg$^{-1}$\n0-%.0f m bulk shear\n\
	    = %.0f kts' % \
	        (plcl.magnitude, SBCAPE.magnitude, alt[-3], bulkshear)
	    boxprops = dict(boxstyle='round', facecolor='none')
	    ax_data.text(0.05, 0.95, datastr, transform=ax_data.transAxes, 
	    	fontsize=14, verticalalignment='top', bbox=boxprops)
	    # Logos
	    ax_png = fignum.add_subplot(gs[3, 3])
	    img = mpimg.imread(logoName)
	    plt.axis('off')
	    plt.imshow(img)
	else:
	    # Logos
	    ax_png = fignum.add_subplot(gs[3, 2:])
	    img = mpimg.imread(logoName)
	    plt.axis('off')
	    plt.imshow(img)

	plt.show(block=False)
	return