# Import all required packages.
from __future__ import division
from astropy import units as u
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import matplotlib.pylab as lab
import numpy as np
from xlrd import open_workbook
import ipdb

# Opening Excel File
wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx') # laptop and work machine
#wb = open_workbook('/home/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx') # home desktop

# Some functions for converting
def ra_to_deg(ra):
    pieces = ra.split(' ')
    hrs = pieces[0]
    mins = pieces[1]
    secs = pieces[2]
    return float(hrs)*(360/24) + float(mins)*(15./60.) + float(secs)*(15./3600.)

def ra_from_deg(ra):
    tmp_hrs = ra/15.
    hrs = np.floor(tmp_hrs)
    tmp_mins = (tmp_hrs-hrs)*60
    mins = np.floor(tmp_mins)
    tmp_secs = (tmp_mins-mins)*60
    secs = round(tmp_secs,2)
    return "{} {} {}".format(int(hrs),str(int(mins)).zfill(2),str(secs).zfill(5))

def dec_to_deg(dec):
    pieces = dec.split(' ')
    #check to see if unicode minus sign is present
    if dec[0] == u'\u2212':
        degs = -1*float(pieces[0][1:])
        mins = -1*float(pieces[1])
        secs = -1*float(pieces[2])
        return degs+mins/60.+secs/3600.
    # if not, proceed normally
    else:
        degs = float(pieces[0])
        mins = float(pieces[1])
        secs = float(pieces[2])
        return degs+mins/60.+secs/3600.

def dec_from_deg(dec):
    if dec < 0:
        degs = np.ceil(dec)
        tmp_mins = (degs-dec)*60
        mins = np.floor(tmp_mins)
        tmp_secs = (tmp_mins-mins)*60
        secs = np.round(tmp_secs,1)
        return "{} {} {}".format(str(int(degs)).zfill(3),str(int(mins)).zfill(2),str(secs).zfill(4))
    if dec > 0:
        degs = np.floor(dec)
        tmp_mins = (dec-degs)*60
        mins = np.floor(tmp_mins)
        tmp_secs = (tmp_mins-mins)*60
        secs = np.round(tmp_secs,1)
        return "{} {} {}".format("+"+str(int(degs)).zfill(2),str(int(mins)).zfill(2),str(secs).zfill(4))

def color_data(colors,i):
    if method_raw[i].value == 'WL':
        colors.append('purple')
    elif method_raw[i].value == 'SL':
        colors.append('red')
    elif method_raw[i].value == 'WL+SL':
        colors.append('black')
    elif method_raw[i].value == 'LOSVD':
        colors.append('orange')
    elif method_raw[i].value == 'X-ray':
        colors.append('green')
    elif method_raw[i].value == 'CM':
        colors.append('blue')
    return colors

def scatter_data(scatters,i):
    if method_raw[i].value == 'WL':
        scatters.append('d')
    elif method_raw[i].value == 'SL':
        scatters.append('s')
    elif method_raw[i].value == 'WL+SL':
        scatters.append('o')
    elif method_raw[i].value == 'LOSVD':
        scatters.append('^')
    elif method_raw[i].value == 'X-ray':
        scatters.append('*')
    elif method_raw[i].value == 'CM':
        scatters.append('x')
    return scatters

    
# Creating Data Structures
ra_raw = []
dec_raw = []
ra_hms = []
dec_dms = []
ra_deg = []
dec_deg = []
method_raw = []
colors = []
scatters = []
# Retrieving RA/DEC
for sheet in wb.sheets():
    if sheet.name == 'Sheet1':
        for i in range(sheet.nrows):
            if i != 0:
                ra_raw.append(sheet.row(i)[19])
                dec_raw.append(sheet.row(i)[20])
                method_raw.append(sheet.row(i)[2])

# Get only good values
badlist = ['TBD','nan','infty']
for i in range(len(ra_raw)):
    if ra_raw[i].value not in badlist and dec_raw[i].value not in badlist:
        ra_hms.append(ra_raw[i].value)
        dec_dms.append(dec_raw[i].value)
        colors = color_data(colors,i)
        scatters = scatter_data(scatters,i)
                    
# Covnert to degrees
for i in range(len(ra_hms)):
    ra_deg.append(ra_to_deg(ra_hms[i]))
    dec_deg.append(dec_to_deg(dec_dms[i]))

# Transform the data into a SkyCoord object.
c = SkyCoord(ra=ra_deg*u.degree, dec=dec_deg*u.degree, frame='icrs')

# Matplotlib needs the coordinates in radians, so we have to convert them.
# Furtermore matplotlib needs the RA coordinate in the range between -pi
# and pi, not 0 and 2pi.
ra_rad = c.ra.radian
dec_rad = c.dec.radian
ra_rad[ra_rad > np.pi] -= 2. * np.pi

# Now plot the data in projection with a grid.
fig = plt.figure()
allsky = lab.subplot(111,projection="rectilinear")
lab.grid(True)
lab.xticks(range(0,375,15))
for label in allsky.axes.xaxis.get_ticklabels()[::2]:
    label.set_visible(False)
lab.yticks(range(-90,105,15))
#for label in allsky.axes.yaxis.get_ticklabels()[::2]:
#    label.set_visible(False)

#'''
from astropy import units as u
from astropy.coordinates import SkyCoord

l_list = np.linspace(-180,180,150) # galactic plane
b_list = np.zeros(150) # galactic plane
ra_mw, dec_mw = ([],[])
for i in range(len(l_list)):
    coord = SkyCoord(l=l_list[i]*u.degree, b=b_list[i]*u.degree, frame='galactic')
    ra_mw.append(coord.icrs.ra.value)
    dec_mw.append(coord.icrs.dec.value)

for i in range(len(l_list)):
    plt.scatter(ra_mw[i], dec_mw[i], marker='.', color='pink')
#'''


print("Plotting {} unique clusters on the sky...".format(len(ra_rad)))

first_xray,first_wl,first_sl,first_wlsl,first_cm,first_losvd = (True,True,True,True,True,True)
for i in range(len(ra_deg)):
    if scatters[i] == '*' and first_xray is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='X-ray')
        first_xray = False
    elif scatters[i] == 'd' and first_wl is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='WL')
        first_wl = False
    elif scatters[i] == 's' and first_sl is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='SL')
        first_sl = False
    elif scatters[i] == 'o' and first_wlsl is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='WL+SL')
        first_wlsl = False
    elif scatters[i] == 'x' and first_cm is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='CM')
        first_cm = False
    elif scatters[i] == '^' and first_losvd is True:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8,label='LOSVD')
        first_losvd = False
    else:
        plt.plot(ra_deg[i], dec_deg[i], '{}'.format(scatters[i]), color=colors[i], markersize=8, alpha=0.75)

plt.ylabel('Dec (deg)',fontsize=18)
plt.xlabel('RA (deg)',fontsize=18)
plt.legend(loc=0,numpoints=1,fontsize=9)
plt.xlim(0,360)
plt.ylim(-90,90)
plt.show()
