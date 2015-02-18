from xlrd import open_workbook

from astropy.io import ascii

import Uncertainties as un
import MConvert_Personal as mc

import ipdb

# Importing original conc/Mass data in delta_vir=500, from text file
fh = ascii.read('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/DE10_1_orig.txt',
                delimiter=',')

# Importing redshifts from excel db (entries have been made)
wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')
sheet_names = wb.sheet_names()
sheet1 = wb.sheet_by_name(sheet_names[0])
clusters = sheet1.col_values(0)
redshift = sheet1.col_values(1)
short_refs = sheet1.col_values(15)
de10_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'DE10.1']
de10_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'DE10.1']
de10_deltavir = [mc.DeltaFinder(0.3,0.7,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'DE10.1']

# Splitting up the data
de10_m500 = []
de10_m500_p = []
de10_m500_m = []
de10_c500 = []
de10_c500_p = []
de10_c500_m = []

fh_clusters = [fh[1][:][j].strip("'") for j in range(len(fh[0][:]))]
for i in range(len(de10_clusters)):
    for j in range(len(fh_clusters)):
        if fh_clusters[j] == de10_clusters[i]:
            de10_c500.append(fh[j][2])
            de10_c500_p.append(fh[j][3])
            de10_c500_m.append(fh[j][4])
            de10_m500.append(fh[j][5])
            de10_m500_p.append(fh[j][6])
            de10_m500_m.append(fh[j][7])


# testing to make sure that data are in the expected format
# a few orders of business first
assert len(de10_clusters) == len(de10_z)
assert len(de10_z) == len(de10_deltavir)
# make sure both concentration and mass measurements are there in all cases
assert len(de10_m500) == len(de10_c500)
assert 'TBD' not in de10_m500
assert 'TBD' not in de10_c500
# make sure both mass uncertainties are there in all cases, and are the same
assert len(de10_m500_p) == len(de10_m500), "M500 uncertainties appear to be missing for some clusters..."
assert len(de10_m500_p) == len(de10_m500_m), "M500 uncertainties do not have the same length..."
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(de10_c500_p) == len(de10_c500), "c500 uncertainties appear to be missing for some clusters..."
assert len(de10_c500_p) == len(de10_c500_m), "c500 uncertainties do not have the same length..."

# converting delta=500 to delta_vir first; then 200 next
# asymmetric uncertainties will be propagated, by conserving the fractional error
delta_orig = 500
c_accuracy = 2 # number of places after the decimal to round to
m_accuracy = 3 # number of places after the decimal to round to

de10_mvir = [round(mc.Mconvert(de10_m500[i],delta_orig,de10_deltavir[i],de10_c500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_mvir_p = [round(de10_mvir[i]*(de10_m500_p[i]/de10_m500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_mvir_m = [round(de10_mvir[i]*(de10_m500_m[i]/de10_m500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_cvir = [round(mc.Cconvert(de10_m500[i],delta_orig,de10_deltavir[i],de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]
de10_cvir_p = [round(de10_cvir[i]*(de10_c500_p[i]/de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]
de10_cvir_m = [round(de10_cvir[i]*(de10_c500_m[i]/de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]

delta_new = 200
de10_m200 = [round(mc.Mconvert(de10_m500[i],delta_orig,delta_new,de10_c500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_m200_p = [round(de10_m200[i]*(de10_m500_p[i]/de10_m500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_m200_m = [round(de10_m200[i]*(de10_m500_m[i]/de10_m500[i]),m_accuracy) for i in range(len(de10_m500))]
de10_c200 = [round(mc.Cconvert(de10_m500[i],delta_orig,delta_new,de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]
de10_c200_p = [round(de10_c200[i]*(de10_c500_p[i]/de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]
de10_c200_m = [round(de10_c200[i]*(de10_c500_m[i]/de10_c500[i]),c_accuracy) for i in range(len(de10_m500))]

for i in range(len(de10_mvir)):
    print "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(i+1,de10_clusters[i],de10_z[i],de10_c200[i],de10_c200_p[i],de10_c200_m[i],de10_m200[i],de10_m200_p[i],de10_m200_m[i],de10_cvir[i],de10_cvir_p[i],de10_cvir_m[i],de10_mvir[i],de10_mvir_p[i],de10_mvir_m[i])



