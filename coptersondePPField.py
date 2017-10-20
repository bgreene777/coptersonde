import matplotlib
matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
import matplotlib.image as mpimg
import matplotlib.gridspec as gridspec
import matplotlib.dates as mpdates
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.basemap import Basemap
import numpy as np
import csv
import math
from datetime import datetime, timedelta
import pytz
import Tkinter
from tkFileDialog import askopenfilename
import scipy.signal as sc
from metpy.plots import Hodograph, SkewT
import metpy.calc as mcalc
from metpy.units import units
import warnings
import cmocean
import geopy.distance
import wxtools

#############
## Version ##
#############
'''
Updated 6 Oct 2017
Brian Greene
University of Oklahoma
For use with raw OU Coptersonde vertical profile data
'''
###############################
## Required pacakges & files ##
###############################
'''
metpy, cmocean, geopy, mesonetgeoinfo.csv
'''
#######################################
## Change the following directories: ##
#######################################
'''
****** os filesep = '/' and windows filesp = '\\' ******
****** this version formatted for mac os ******
'''
# nextcloud directory
myNextcloud = '/Users/briangreene/Nextcloud/thermo/'
# location of raw .csv and .pos files
dataDirName = myNextcloud + 'data/PBL Transition/'
# location of mesonet location csv
mesocsv = myNextcloud + 'documentation/Mesonet/geoinfo.csv'
# location of Logos.png
logoName = myNextcloud + 'documentation/Logos.png'
# location to save csv output
folderSaveFile = '/Users/briangreene/Documents/Coptersonde/KAEFS/Data/'
# location to save png output
folderSavePNG = '/Users/briangreene/Documents/Coptersonde/KAEFS/Figures/'

##################################
## Setup and Raw File Selection ##
##################################

# Ignore dumb warnings
warnings.filterwarnings("ignore",".*GUI is implemented.*")
warnings.filterwarnings("ignore",".*mean of empty slice.*")
warnings.filterwarnings("ignore",".*invalid value encountered in less.*")

# Select file
autofile = raw_input('>>>Use latest CSV? y/n ')
while autofile != 'y' and autofile != 'n':
    autofile = raw_input('>>>Use latest CSV? y/n ')

if autofile == 'y':
    filename = wxtools.findLatestCSV(dataDirName)
elif autofile == 'n':
    autodir = raw_input('>>>Use latest directory? y/n ')
    while autodir != 'y' and autodir != 'n':
        autodir = raw_input('>>>Use latest directory? y/n ')
    if autodir == 'y':
        PBLdir = wxtools.findLatestDir(dataDirName)
        print '>>>Please select CSV.'
        root = Tkinter.Tk()
        root.withdraw()
        root.update()
        filename = askopenfilename(initialdir=PBLdir)
        root.destroy()
    else:
        print '>>>Please select CSV.'
        root = Tkinter.Tk()
        root.withdraw()
        root.update()
        filename = askopenfilename(initialdir = dataDirName)
        root.destroy()
else:
    print 'No file selected!'

fname = filename.rsplit('/', 1)[-1]
print 'CSV file selected: %s' % fname

# Copter Number
copterNum = fname[11]
#copterNum = '1'

# Ask if ppk
dgpsFlag = raw_input('>>>PPK? y/n ')
while dgpsFlag != 'y' and dgpsFlag != 'n':
    dgpsFlag = raw_input('>>>PPK? y/n ')

if dgpsFlag == 'y':
    dgpsFlag = 1
    ppkname = wxtools.findDGPSfile(filename)
    print 'PPK file selected: %s' % ppkname.split('/')[-1]
elif dgpsFlag == 'n':
    dgpsFlag = 0
    print 'No PPK. '
else:
    print '>>>ERROR '

######################
## Import DGPS Data ##
######################

if dgpsFlag:
    dgpsdata = wxtools.dgpsread(ppkname)
    date_dgps = dgpsdata[0].tolist()
    time_dgps = dgpsdata[1].tolist()
    lat_dgps = dgpsdata[2].astype(np.float)
    lon_dgps = dgpsdata[3].astype(np.float)
    alt_dgps = dgpsdata[4].astype(np.float)
    alt_dgps_agl = dgpsdata[5].astype(np.float)

    # Convert time to datetime object
    datestr_dgps = [x+' '+y for x,y in zip(date_dgps,time_dgps)]
    datetime_dgps = [datetime.strptime(idatestr, 
        '%Y/%m/%d %H:%M:%S.%f') for idatestr in datestr_dgps]
    for i in range(len(datetime_dgps)):
        datetime_dgps[i] = datetime_dgps[i].replace(tzinfo=pytz.utc)

#################
## Flight data ##
#################

deltaZ = 10.
# Automatically find max height
if dgpsFlag:
    maxHeight = round(np.nanmax(alt_dgps_agl), -1) + deltaZ
    print 'Maximum height AGL (PPK): %5.4f m' % np.nanmax(alt_dgps_agl)
else:
    H = []
    a = open(filename)
    reader = csv.DictReader(a)
    for line in reader:
        H.append(float(line['Altitude']))

    a.close()
    print 'Maximum height AGL (copter): %5.4f m' % np.nanmax(H)
    maxHeight = round(np.nanmax(H), -1) + deltaZ

nHeights = int(maxHeight / deltaZ)
sampleHeights_m = np.linspace(10., maxHeight, num=nHeights)

#################
## Import data ##
#################

copterdata = wxtools.csvread_raw(filename)
time = copterdata[0]
lat = copterdata[1]
lon = copterdata[2]
alt = copterdata[3]
p = copterdata[4]
roll = copterdata[5]
pitch = copterdata[6]
yaw = copterdata[7]
RH1 = copterdata[8]
RH2 = copterdata[9]
RH3 = copterdata[10]
RH4 = copterdata[11]
RHT1 = copterdata[12]
RHT2 = copterdata[13]
RHT3 = copterdata[14]
RHT4 = copterdata[15]
T1 = copterdata[16]
T2 = copterdata[17]
T3 = copterdata[18]
T4 = copterdata[19]
gspd = copterdata[20]
vspd = copterdata[21]
ax = copterdata[22]
ay = copterdata[23]

# Convert time to datetime object
dt_copter = [datetime.utcfromtimestamp(i) for i in time]
# Apply UTC timezone
for i in range(len(dt_copter)):
    dt_copter[i] = dt_copter[i].replace(tzinfo=pytz.utc)

###################
## Site location ##
###################

# isfindsite = raw_input('>>>Automatically find Mesonet site? y/n ')
# while isfindsite != 'y' and isfindsite != 'n':
#     isfindsite = raw_input('>>>Automatically find Mesonet site? y/n ')
isfindsite = 'y'

if isfindsite == 'y':
    meso = wxtools.findSite(lat[0], lon[0])
    sitename = meso.split('/')[0]
    sitelong = meso.split('/')[1]
elif isfindsite == 'n':
    sitename = raw_input('>>>Please enter 4-digit site name: ')
    sitelong = sitename
else:
    print '>>>ERROR'

################################################
## Apply offsets - will do calibrations later ##
################################################

#################################################
## Apply Median Filter - preserves same length ##
#################################################

T1 = sc.medfilt(T1, 3)
T2 = sc.medfilt(T2, 3)
T3 = sc.medfilt(T3, 3)
T4 = sc.medfilt(T4, 3)
RH1 = sc.medfilt(RH1, 3)
RH2 = sc.medfilt(RH2, 3)
RH3 = sc.medfilt(RH3, 3)
RH4 = sc.medfilt(RH4, 3)

#####################
## Calculate Winds ##
#####################

nVals = len(roll)

psi_deg = np.zeros(nVals)
az_deg = np.zeros(nVals)
for j in range(nVals):
    crol = np.cos(roll[j] * np.pi / 180.)
    srol = np.sin(roll[j] * np.pi / 180.)
    cpit = np.cos(pitch[j] * np.pi / 180.)
    spit = np.sin(pitch[j] * np.pi / 180.)
    cyaw = np.cos(yaw[j] * np.pi / 180.)
    syaw = np.sin(yaw[j] * np.pi / 180.)

    Rx = np.matrix( [[1,0,0], [0,crol,srol], [0,-srol,crol]] )
    Ry = np.matrix( [[cpit,0,-spit], [0,1,0], [spit,0,cpit]] )
    Rz = np.matrix( [[cyaw,-syaw,0], [syaw,cyaw,0], [0,0,1]] )
    R = Rz * Ry * Rx

    psi_deg[j] = np.arccos(R[2,2]) * 180./np.pi
    az_deg[j] = np.arctan2(R[1,2],R[0,2]) * 180./np.pi

Speed_mps = 34.5 * np.sqrt(np.tan(psi_deg * np.pi/180.)) - 6.4
Direction_deg = az_deg
Speed_mps[Speed_mps<0.] = np.nan

# Fix negative values
iNeg = np.squeeze(np.where(Direction_deg < 0.))
Direction_deg[iNeg] = Direction_deg[iNeg] + 360.

###########################################
## Automatically Determine Starting Time ##
###########################################

updown = raw_input('>>>Enter 1 for ascent or 2 for descent: ')
if updown == '1':
    profup = 1
    profdown = 0
    print 'Selected: ascent'
elif updown == '2':
    profup = 0
    profdown = 1
    print 'Selected: descent'
else:
    print 'uhhhhh'


t_copter = mpdates.date2num(dt_copter)

if profup:
    timeTakeoff = wxtools.findStart(vspd, alt, t_copter)
    if dgpsFlag:
        timeMax = datetime_dgps[np.argmax(alt_dgps)]
        t_dgps = mpdates.date2num(datetime_dgps)

        # Find time diff between copter and dgps
        timeMaxCopter = dt_copter[np.argmax(alt)]
        timeDiff = timeMax - timeMaxCopter
        print 'Time offset: %5.2f seconds' % timeDiff.total_seconds()

        # Correct copter time to match dgps
        dt_copter = [time + timeDiff for time in dt_copter]
        t_copter = mpdates.date2num(dt_copter)
        timeTakeoff += timeDiff

        fig1, ax1 = plt.subplots(1)
        plt.plot(t_dgps,alt_dgps)
        plt.title('Copter & DGPS Altitude - Takeoff & Max in Red')
        plt.ylabel('Copter & DGPS Altitude (AGL, ASL) (m)')

    else:
        timeMax = dt_copter[np.argmax(alt)]
        fig1, ax1 = plt.subplots(1)
        plt.title('Altitude vs. Time - Takeoff and Max Times Indicated in Red')
        plt.ylabel('Copter Altitude AGL (m)')

    plt.plot(t_copter, alt)
    plt.plot(timeTakeoff, 10, 'r.')
    plt.plot(timeMax, np.max(alt), 'r.')
    plt.xlabel('Time UTC - Copter')
    ax1.xaxis.set_major_locator(mpdates.MinuteLocator(interval=5))
    ax1.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M:%S'))

    print 'Takeoff Time (UTC): %s' % timeTakeoff.isoformat(' ')
    print 'Max Time (UTC): %s' % timeMax.isoformat(' ')

    plt.show(block=False)

if profdown:
    # "timeTakeoff" = time of apex of flight for descending profile
    if dgpsFlag:
        timeTakeoff = datetime_dgps[np.argmax(alt_dgps)]
        t_dgps = mpdates.date2num(datetime_dgps)

        # Find time diff between copter and dgps
        timeMaxCopter = dt_copter[np.argmax(alt)]
        timeDiff = timeTakeoff - timeMaxCopter
        print 'Time offset: %5.2f seconds' % timeDiff.total_seconds()

        # Correct copter time to match dgps
        dt_copter = [time + timeDiff for time in dt_copter]
        t_copter = mpdates.date2num(dt_copter)
        timeEnd = wxtools.findEnd(vspd, alt, t_copter)

        fig1, ax1 = plt.subplots(1)
        plt.plot(t_dgps,alt_dgps)
        plt.title('Copter & DGPS Altitude - Takeoff & Max in Red')
        plt.ylabel('Copter & DGPS Altitude (AGL, ASL) (m)')

    else:
        timeTakeoff = dt_copter[np.argmax(alt)]
        timeEnd = wxtools.findEnd(vspd, alt, t_copter)
        fig1, ax1 = plt.subplots(1)
        plt.title('Altitude vs. Time - Takeoff and Max Times Indicated in Red')
        plt.ylabel('Copter Altitude AGL (m)')

    plt.plot(t_copter, alt)
    plt.plot(timeTakeoff, np.max(alt), 'r.')
    plt.plot(timeEnd, 10, 'r.')
    plt.xlabel('Time UTC - Copter')
    ax1.xaxis.set_major_locator(mpdates.MinuteLocator(interval=5))
    ax1.xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M:%S'))

    print 'Takeoff Time (UTC): %s' % timeTakeoff.isoformat(' ')
    print 'End Time (UTC): %s' % timeEnd.isoformat(' ')

    plt.show(block=False)


#################################################
## Find indices for coptersonde and dgps files ##
#################################################

indCopter = []
count = 0
if profup:
    for d in dt_copter:
        if (d >= timeTakeoff) & (d <= timeMax):
            indCopter.append(count)
        count += 1
    indCopter = np.array(indCopter)

    if dgpsFlag:
        indDGPS = []
        count = 0
        for d in datetime_dgps:
            if (d >= timeTakeoff) & (d <= timeMax):
                indDGPS.append(count)
            count += 1
        indDGPS = np.array(indDGPS)
if profdown:
    for d in dt_copter:
        if (d >= timeTakeoff) & (d <= timeEnd):
            indCopter.append(count)
        count += 1
    indCopter = np.array(indCopter)

    if dgpsFlag:
        indDGPS = []
        count = 0
        for d in datetime_dgps:
            if (d >= timeTakeoff) & (d <= timeEnd):
                indDGPS.append(count)
            count += 1
        indDGPS = np.array(indDGPS)

#########################################
## Plot filtered data in chosen Domain ##
#########################################

fig4, axarr = plt.subplots(2, sharex=True)
axarr[0].plot(t_copter[indCopter], T1[indCopter], label='T1')
axarr[0].plot(t_copter[indCopter], T2[indCopter], label='T2')
axarr[0].plot(t_copter[indCopter], T3[indCopter], label='T3')
axarr[0].plot(t_copter[indCopter], T4[indCopter], label='T4')
axarr[0].set_title('Temperature Median Filtered')
axarr[0].xaxis.set_major_locator(mpdates.MinuteLocator(interval=1))
axarr[0].xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
axarr[0].legend(loc=0)

axarr[1].plot(t_copter[indCopter], RH1[indCopter], label='RH1')
axarr[1].plot(t_copter[indCopter], RH2[indCopter], label='RH2')
axarr[1].plot(t_copter[indCopter], RH3[indCopter], label='RH3')
axarr[1].plot(t_copter[indCopter], RH4[indCopter], label='RH4')
axarr[1].set_title('RH Median Filtered')
axarr[1].xaxis.set_major_locator(mpdates.MinuteLocator(interval=1))
axarr[1].xaxis.set_major_formatter(mpdates.DateFormatter('%H:%M'))
axarr[1].legend(loc=0)

plt.show(block=False)

###################################
## Prompt to select data to omit ##
###################################

cT = input('>>>Enter T sensor #s separated by commas to omit,\
 0 if none, or 5 if all: ')
cH = input('>>>Enter RH sensor #s separated by commas to omit,\
 0 if none, or 5 if all: ')

if not isinstance(cT, list):
    cT = np.array(cT)
    np.append(cT,0)
if not isinstance(cH, list):
    cH = np.array(cH)
    np.append(cH,0)

if 1 in cT:
    T1 = np.array([np.nan for i in T1])
if 2 in cT:
    T2 = np.array([np.nan for i in T2])
if 3 in cT:
    T3 = np.array([np.nan for i in T3])
if 4 in cT:
    T4 = np.array([np.nan for i in T4])
if 5 in cT:
    T1, T2, T3, T4 = (np.array([np.nan for i in T1]) for j in range(4))
if 1 in cH:
    RH1 = np.array([np.nan for i in RH1])
if 2 in cH:
    RH2 = np.array([np.nan for i in RH2])
if 3 in cH:
    RH3 = np.array([np.nan for i in RH3])
if 4 in cH:
    RH4 = np.array([np.nan for i in RH4])
if 5 in cH:
    RH1, RH2, RH3, RH4 = (np.array([np.nan for i in RH1]) for j in range(4))

#####################################
## Average the values over delta z ##
#####################################

t1, t2, t3, t4, rh1, rh2, rh3, rh4, pres, wind, direction = (np.array(
    [np.nan for i in sampleHeights_m]) for j in range(11))

if dgpsFlag:
    datetime_dgps = np.array(datetime_dgps)

dt_copter = np.array(dt_copter)
for iHeight in np.arange(nHeights):
    if dgpsFlag:
        ind = np.squeeze(np.where(
            (sampleHeights_m[iHeight] - deltaZ/2. <= alt_dgps_agl[indDGPS])\
            & (alt_dgps_agl[indDGPS] <= sampleHeights_m[iHeight] + deltaZ/2.)))
        if ind.size != 0:
            ihere = 1
            indTimeAvg = []
            count = 0
            for d in dt_copter[indCopter]:
                if (d >= datetime_dgps[indDGPS][ind][0]) &\
                (d <= datetime_dgps[indDGPS][ind][-1]):
                    indTimeAvg.append(count)
                count += 1
        else:
            ihere = 0

    else:
        indTimeAvg = np.squeeze(np.where(
                (sampleHeights_m[iHeight] - deltaZ/2. <= alt[indCopter]) &\
                (alt[indCopter] <= sampleHeights_m[iHeight] + deltaZ/2.)) )
        ihere = 1

    indTimeAvg = np.array(indTimeAvg)

    if ihere:
        t1[iHeight] = np.nanmean(T1[indCopter][indTimeAvg], 0)
        t2[iHeight] = np.nanmean(T2[indCopter][indTimeAvg], 0)
        t3[iHeight] = np.nanmean(T3[indCopter][indTimeAvg], 0)
        t4[iHeight] = np.nanmean(T4[indCopter][indTimeAvg], 0)

        rh1[iHeight] = np.nanmean(RH1[indCopter][indTimeAvg], 0)
        rh2[iHeight] = np.nanmean(RH2[indCopter][indTimeAvg], 0)
        rh3[iHeight] = np.nanmean(RH3[indCopter][indTimeAvg], 0)
        rh4[iHeight] = np.nanmean(RH4[indCopter][indTimeAvg], 0)

        pres[iHeight] = np.nanmean(p[indCopter][indTimeAvg], 0)

        wind[iHeight] = np.nanmean(Speed_mps[indCopter][indTimeAvg], 0)
        direction[iHeight] = np.arctan2(
            np.nanmean(np.sin(Direction_deg[indCopter][indTimeAvg]
                * np.pi / 180.), 0),
            np.nanmean(np.cos(Direction_deg[indCopter][indTimeAvg]
                * np.pi / 180.), 0) ) * 180. / np.pi

########################################
## Convert wind to kts, fix direction ##
########################################

#wind_kts = np.array([w * 1.94 for w in wind])
wind_kts = wind * 1.94

iNeg = np.squeeze(np.where(direction < 0.))
direction[iNeg] = direction[iNeg] + 360.

u,v = mcalc.get_wind_components(wind_kts * units.kts, direction * units.deg)
u = u.to(units.kts)
v = v.to(units.kts)

##################################
## Find average between sensors ##
##################################

TArr = np.array([t1, t2, t3, t4])
RHArr = np.array([rh1, rh2, rh3, rh4])

Tmean = np.nanmean(TArr, 0)
RHmean = np.nanmean(RHArr, 0)
if np.isnan(RHmean).all():
    isRH = 0
else:
    isRH = 1

######################################################
## Calculate thermodymanic and convective variables ##
######################################################

theta = np.array(mcalc.potential_temperature(pres * units.mbar,
    (Tmean + 273.15) * units.kelvin))
Td = np.array(mcalc.dewpoint_rh(Tmean * units.degC, RHmean / 100.))
ws = np.array(mcalc.saturation_mixing_ratio(pres * units.mbar,
    (Tmean + 273.15) * units.kelvin))
w = np.multiply(RHmean / 100., ws * 1000.)

# check for RH
if isRH:
    lclpres, lcltemp = mcalc.lcl(pres[0] * units.mbar, 
        Tmean[0] * units.degC, Td[0] * units.degC)
    print 'LCL Pressure: %5.2f %s' % (lclpres.magnitude, lclpres.units)
    print 'LCL Temperature: %5.2f %s' % (lcltemp.magnitude, lcltemp.units)

    # parcel profile
    # determine if there are points sampled above lcl
    ilcl = np.squeeze(np.where((pres * units.mbar) <= lclpres))
    # if not, entire profile dry adiabatic
    if ilcl.size == 0:
        prof = mcalc.dry_lapse(pres * units.mbar, 
            Tmean[0] * units.degC).to('degC')
        isbelowlcl = 1
    # if there are, need to concat dry & moist profile ascents
    else:
        ilcl = ilcl[0]
        prof_dry = mcalc.dry_lapse(pres[:ilcl] * units.mbar,
            Tmean[0] * units.degC).to('degC')
        prof_moist = mcalc.moist_lapse(pres[ilcl:] * units.mbar,
            prof_dry[-1]).to('degC')
        prof = np.concatenate((prof_dry, prof_moist)) * units.degC
        isbelowlcl = 0

    # CAPE
    SBCAPE = wxtools.uavCAPE(Tmean * units.degC, prof, pres * units.hPa)
    print 'Estimated SBCAPE: %4.2f %s' % (SBCAPE.magnitude, SBCAPE.units)
else:
    isbelowlcl = 0

# Wind shear
bulkshear = wind_kts[-3] - wind_kts[0]
print '0-%d m Bulk Shear: %.0f kts' % (sampleHeights_m[-3], bulkshear)

#########################
## Import Mesonet Data ##
#########################
mesoData = wxtools.getMesoData(timeTakeoff.year, timeTakeoff.month,
    timeTakeoff.day, sitename)
if mesoData.size == 0:
    print 'No internet connection detected.'
    RHmeso = np.nan
    T2meso = np.nan
    T9meso = np.nan
    umeso = np.nan
    vmeso = np.nan
    pmeso = np.nan
    Td2meso = np.nan
else:
    print 'Internet connection successful! Pulling Mesonet data...'
    iMesoTime = wxtools.findClosestMesoTime(timeTakeoff)
    tmeso = mesoData[0, iMesoTime]
    tmeso = datetime(timeTakeoff.year, timeTakeoff.month, timeTakeoff.day) + \
        timedelta(minutes=tmeso)

    RHmeso = mesoData[1, iMesoTime]
    T2meso = mesoData[2, iMesoTime]
    T9meso = mesoData[3, iMesoTime]
    umeso = mesoData[4, iMesoTime]
    vmeso = mesoData[5, iMesoTime]
    pmeso = mesoData[6, iMesoTime]
    Td2meso = np.array(mcalc.dewpoint_rh(T2meso * units.degC, RHmeso / 100.))

######################
## Create SkewTLogP ##
######################

print 'Plotting...'
fig5 = plt.figure(figsize=(12,9))
gs = gridspec.GridSpec(4, 4)
skew = SkewT(fig5, rotation=20, subplot=gs[:, :2])

skew.plot(pres, Tmean, 'r', linewidth = 2)
skew.plot(pres, Td, 'g', linewidth = 2)
skew.plot_barbs(pres[0::4], u[0::4], v[0::4], x_clip_radius = 0.12, \
    y_clip_radius = 0.12)

# Plot mesonet surface data and winds
skew.plot(pmeso, T2meso, 'k*', linewidth=2, label='Mesonet 2m T')
skew.plot(pres[0], T9meso, 'r*', linewidth=2, label='Mesonet 9m T')
skew.plot(pmeso, Td2meso, 'g*', linewidth=2, label='Mesonet 2m Td')
skew.plot_barbs(pmeso, umeso, vmeso, barbcolor='r')

plt.legend(loc=4)

# Plot convective parameters
if isRH:
    skew.plot(lclpres, lcltemp, 'ko', markerfacecolor='black')
    skew.plot(pres, prof, 'k', linewidth=2)

# set up plot limits and labels - use LCL as max if higher than profile
if isRH:
    xmin = math.floor(np.nanmin(Td))
else:
    xmin = math.floor(np.nanmin(Tmean))
xmax = math.floor(np.nanmax(Tmean)) + 20
if isbelowlcl:
    ymin = round((lclpres / units.mbar), -1) - 10
else:
    ymin = round(np.nanmin(pres),-1) - 10
    
ymax = round(np.nanmax(pres),-1) + 10

skew.ax.set_ylim(ymax,ymin)
skew.ax.set_xlim(xmin,xmax)
skew.ax.set_yticks(np.arange(ymin,ymax+10,10))

skew.ax.set_xlabel('Temperature ($^\circ$C)')
skew.ax.set_ylabel('Pressure (hPa)')
titleName = 'Coptersonde-%s %s UTC - %s' % (copterNum, 
    timeTakeoff.strftime('%d-%b-%Y %H:%M:%S'), sitename)
skew.ax.set_title(titleName)

skew.plot_dry_adiabats(linewidth=0.75)
skew.plot_moist_adiabats(linewidth=0.75)
skew.plot_mixing_lines(linewidth=0.75)

# Hodograph
ax_hod = fig5.add_subplot(gs[:2,2:])
#gs.tight_layout(fig5)
if np.nanmax(wind_kts) > 18:
    comprange = 35
else:
    comprange = 20

h = Hodograph(ax_hod, component_range=comprange)
h.add_grid(increment=5)
h.plot_colormapped(u,v,pres, cmap=cmocean.cm.deep_r)
ax_hod.set_title('Hodograph (kts)')
ax_hod.yaxis.set_ticklabels([])
#ax_hod.set_xlabel('Wind Speed (kts)')

# Map - Oklahoma
llcrnrlat = 33.6
urcrnrlat = 37.2
llcrnrlon = -103.2
urcrnrlon = -94.2
ax_map = fig5.add_subplot(gs[2, 2:])

m = Basemap(projection='merc', llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat, 
    llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon, lat_ts=20, resolution='l',
    ax=ax_map)

print 'Basemap...'
m.drawcounties()
m.drawstates()
x,y = m(lon[0], lat[0])
plt.plot(x,y,'b.')
plt.text(x+40000, y-5000, sitelong, bbox=dict(facecolor='yellow', alpha=0.5))

if isRH:
    # Convective parameter values
    ax_data = fig5.add_subplot(gs[3, 2])
    plt.axis('off')
    datastr = 'LCL = %.0f hPa\nFake News CAPE = %.0f J kg$^{-1}$\n0-%.0f m bulk shear\n\
        = %.0f kts' % \
        (lclpres.magnitude, SBCAPE.magnitude, sampleHeights_m[-3], bulkshear)
    boxprops = dict(boxstyle='round', facecolor='none')
    ax_data.text(0.05, 0.95, datastr, transform=ax_data.transAxes, fontsize=14,
        verticalalignment='top', bbox=boxprops)
    # Logos
    ax_png = fig5.add_subplot(gs[3, 3])
    img = mpimg.imread(logoName)
    plt.axis('off')
    plt.imshow(img)
else:
    # Logos
    ax_png = fig5.add_subplot(gs[3, 2:])
    img = mpimg.imread(logoName)
    plt.axis('off')
    plt.imshow(img)

plt.show(block=False)

######################
## Save csv and png ##
######################

s = raw_input('>>>Save csv and figures? y/n ')
while s != 'y' and s != 'n':
    s = raw_input('>>>Save csv and figures? y/n ')
if s == 'y':
    saveFileName = folderSaveFile + '%s/%s-OUUAS%s-%s.csv' % \
        (timeTakeoff.strftime('%Y%m%d'), timeTakeoff.strftime('%Y%m%d_%H%M%S'),
            copterNum, sitename)
    headers = ('Lat', 'Lon', 'AltAGL(m)', 'p(hPa)', 'T(C)', 'Td(C)',\
        'RH(percent)', 'w(gKg-1)', 'Theta(K)', 'Speed(ms-1)', 'Dir(deg)')
    fw = open(saveFileName,'wb')
    writer = csv.writer(fw,delimiter=',')
    writer.writerow(headers)
    for i in range(nHeights):
        if (i == nHeights-1) & np.isnan(Tmean[i]):
            print 'Reached end of csv'
            break
        else:
            writer.writerow( (lat[0], lon[0], sampleHeights_m[i],
                round(pres[i], 2), round(Tmean[i], 2), round(Td[i], 2),
                round(RHmean[i], 2), round(w[i], 2), round(theta[i], 2),
                round(wind[i], 2), round(direction[i], 2)) )

    fw.close()
    print 'Finished saving %s' % saveFileName.split('/')[-1]

    # png
    saveFileNamePNG = folderSavePNG + '%s/%s-OUUAS%s-%s.png' % \
        (timeTakeoff.strftime('%Y%m%d'), timeTakeoff.strftime('%Y%m%d_%H%M%S'),
            copterNum, sitename)
    fig5.savefig(saveFileNamePNG)
    print 'Finished saving %s' % saveFileNamePNG.split('/')[-1]

elif s == 'n':
    print 'Files not saved. '
else:
    print '>>>How did you get here??'

plt.show(block=False)

#####################
## Quit when ready ##
#####################

q = raw_input('>>>Press enter to quit. ')
while q != '':
    q = raw_input('>>>Press enter to quit. ')

plt.close('all')
print 'Post Processing Complete.'
print '               ___________                   '
print '              /    O      \\                 '
print '<====>       /      U      \\       <====>   '
print '  []________/_______________\\________[]     '
print '           ||---------------||               '
print '           ||               ||               '
print '           ||               ||               '

