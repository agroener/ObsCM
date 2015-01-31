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
wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')

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

# Creating Data Structures
ra_raw = []
dec_raw = []
ra_hms = []
dec_dms = []
ra_deg = []
dec_deg = []

# Retrieving RA/DEC
for sheet in wb.sheets():
    if sheet.name == 'Sheet1':
        for i in range(sheet.nrows):
            if i != 0:
                ra_raw.append(sheet.row(i)[19])
                dec_raw.append(sheet.row(i)[20])

# Get only good values
badlist = ['TBD','nan','infty']
for i in range(len(ra_raw)):
    if ra_raw[i].value not in badlist and dec_raw[i].value not in badlist:
        ra_hms.append(ra_raw[i].value)
        dec_dms.append(dec_raw[i].value)            
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

# Now plot the data in Aitoff projection with a grid.
#fig = plt.figure()
#lab.subplot(111,projection="aitoff")
#lab.title(r"Aitoff Projection")
#lab.grid(True)
#plt.plot(ra_rad, dec_rad, 'o')
#plt.show()
