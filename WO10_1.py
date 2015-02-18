from xlrd import open_workbook

from astropy.io import ascii

import Uncertainties as un
import MConvert_Personal as mc

import ipdb

# Importing original conc/Mass data in delta_vir=102, from text file
fh = ascii.read('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/WO10_1_orig.txt',
                delimiter=',')

# Importing redshifts from excel db (entries have been made)
wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')
sheet_names = wb.sheet_names()
sheet1 = wb.sheet_by_name(sheet_names[0])
clusters = sheet1.col_values(0)
redshift = sheet1.col_values(1)
short_refs = sheet1.col_values(15)
wo10_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'WO10.1']
wo10_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'WO10.1']
wo10_deltavir = [mc.DeltaFinder(0.3,0.7,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'WO10.1']

# Splitting up the data
wo10_m102 = []
wo10_m102_p = []
wo10_m102_m = []
wo10_c102 = []
wo10_c102_p = []
wo10_c102_m = []

fh_clusters = [fh[0][:][j].strip("'") for j in range(len(fh[0][:]))]

for i in range(len(wo10_clusters)):
    for j in range(len(fh_clusters)):
        if fh_clusters[j] == wo10_clusters[i]:
            wo10_m102.append(fh[j][1])
            wo10_m102_p.append(fh[j][2])
            wo10_m102_m.append(fh[j][3])
            wo10_c102.append(fh[j][4])
            wo10_c102_p.append(fh[j][5])
            wo10_c102_m.append(fh[j][6])

# testing to make sure that data are in the expected format
# a few orders of business first
assert len(wo10_clusters) == len(wo10_z)
assert len(wo10_z) == len(wo10_deltavir)
# make sure both concentration and mass measurements are there in all cases
assert len(wo10_m102) == len(wo10_c102)
assert 'TBD' not in wo10_m102
assert 'TBD' not in wo10_c102
# make sure both mass uncertainties are there in all cases, and are the same
assert len(wo10_m102_p) == len(wo10_m102), "M102 uncertainties appear to be missing for some clusters..."
assert len(wo10_m102_p) == len(wo10_m102_m), "M102 uncertainties do not have the same length..."
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(wo10_c102_p) == len(wo10_c102), "c102 uncertainties appear to be missing for some clusters..."
assert len(wo10_c102_p) == len(wo10_c102_m), "c102 uncertainties do not have the same length..."

# converting delta=102 to delta_vir first; then 200 next
# asymmetric uncertainties will be propagated, by conserving the fractional error
delta_orig = 102
accuracy = 2 # number of places after the decimal to round to

wo10_mvir = [round(mc.Mconvert(wo10_m102[i],delta_orig,wo10_deltavir[i],wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_mvir_p = [round(wo10_mvir[i]*(wo10_m102_p[i]/wo10_m102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_mvir_m = [round(wo10_mvir[i]*(wo10_m102_m[i]/wo10_m102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_cvir = [round(mc.Cconvert(wo10_m102[i],delta_orig,wo10_deltavir[i],wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_cvir_p = [round(wo10_cvir[i]*(wo10_c102_p[i]/wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_cvir_m = [round(wo10_cvir[i]*(wo10_c102_m[i]/wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]

delta_new = 200
wo10_m200 = [round(mc.Mconvert(wo10_m102[i],delta_orig,delta_new,wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_m200_p = [round(wo10_m200[i]*(wo10_m102_p[i]/wo10_m102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_m200_m = [round(wo10_m200[i]*(wo10_m102_m[i]/wo10_m102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_c200 = [round(mc.Cconvert(wo10_m102[i],delta_orig,delta_new,wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_c200_p = [round(wo10_c200[i]*(wo10_c102_p[i]/wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]
wo10_c200_m = [round(wo10_c200[i]*(wo10_c102_m[i]/wo10_c102[i]),accuracy) for i in range(len(wo10_m102))]


for i in range(len(wo10_mvir)):
    print "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(i+1,wo10_clusters[i],wo10_z[i],wo10_c200[i],wo10_c200_p[i],wo10_c200_m[i],wo10_m200[i],wo10_m200_p[i],wo10_m200_m[i],wo10_cvir[i],wo10_cvir_p[i],wo10_cvir_m[i],wo10_mvir[i],wo10_mvir_p[i],wo10_mvir_m[i])
