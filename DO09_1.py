from xlrd import open_workbook

from astropy.table import Table
from astropy.io import ascii

import Uncertainties as un
import MConvert_Personal as mc

import ipdb

# Opening Excel File
wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')

# First sheet
sheet_names = wb.sheet_names()
sheet1 = wb.sheet_by_name(sheet_names[0])

# Get headers and data
headers = []

clusters = sheet1.col_values(0)
headers.append(clusters.pop(0))

redshift = sheet1.col_values(1)
headers.append(redshift.pop(0))

methods = sheet1.col_values(2)
headers.append(methods.pop(0))

c200 = sheet1.col_values(3)
headers.append(c200.pop(0))

c200_plus = sheet1.col_values(4)
headers.append(c200_plus.pop(0))

c200_minus = sheet1.col_values(5)
headers.append(c200_minus.pop(0))

m200 = sheet1.col_values(6)
headers.append(m200.pop(0))

m200_plus = sheet1.col_values(7)
headers.append(m200_plus.pop(0))

m200_minus = sheet1.col_values(8)
headers.append(m200_minus.pop(0))

cvir = sheet1.col_values(9)
headers.append(cvir.pop(0))

cvir_plus = sheet1.col_values(10)
headers.append(cvir_plus.pop(0))

cvir_minus = sheet1.col_values(11)
headers.append(cvir_minus.pop(0))

mvir = sheet1.col_values(12)
headers.append(mvir.pop(0))

mvir_plus = sheet1.col_values(13)
headers.append(mvir_plus.pop(0))

mvir_minus = sheet1.col_values(14)
headers.append(mvir_minus.pop(0))

short_refs = sheet1.col_values(15)
headers.append(short_refs.pop(0))

orig_convention = sheet1.col_values(16)
headers.append(orig_convention.pop(0))

cosmology = sheet1.col_values(17)
headers.append(cosmology.pop(0))

# collect relavent data; there are no nans, TBD's, 
# but there are asymmetric uncertainties here to deal with
do09_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'DO09.1']
do09_m200 = [m200[i] for i in range(len(m200)) if short_refs[i] == 'DO09.1']
do09_m200_p = [m200_plus[i] for i in range(len(m200_plus)) if short_refs[i] == 'DO09.1']
do09_m200_m = [m200_minus[i] for i in range(len(m200_minus)) if short_refs[i] == 'DO09.1']
do09_c200 = [c200[i] for i in range(len(c200)) if short_refs[i] == 'DO09.1']
do09_c200_p = [c200_plus[i] for i in range(len(c200_plus)) if short_refs[i] == 'DO09.1']
do09_c200_m = [c200_minus[i] for i in range(len(c200_minus)) if short_refs[i] == 'DO09.1']
do09_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'DO09.1']
do09_deltavir = [mc.DeltaFinder(0.3,0.7,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'DO09.1']

# testing to make sure that data are in the expected format
# a few orders of business first
assert len(do09_clusters) == len(do09_z)
assert len(do09_z) == len(do09_deltavir)
# make sure both concentration and mass measurements are there in all cases
assert len(do09_m200) == len(do09_c200)
assert 'TBD' not in do09_m200
assert 'TBD' not in do09_c200
# make sure both mass uncertainties are there in all cases, though they are not the same
assert len(do09_m200_p) == len(do09_m200), "M200 uncertainties appear to be missing for some clusters..."
assert len(do09_m200_p) == len(do09_m200_m), "M200 uncertainties do not have the same length..."
# make sure both concentration uncertainties are there in all cases, though they are not the same
assert len(do09_c200_p) == len(do09_c200), "c200 uncertainties appear to be missing for some clusters..."
assert len(do09_c200_p) == len(do09_c200_m), "c200 uncertainties do not have the same length..."

# time to convert; for asymmetric uncertainties, I've decided to keep the fractional error constant
delta_old = 200 # constant for all clusters
accuracy = 2 # number of places after the decimal to round to
do09_mvir = [round(mc.Mconvert(do09_m200[i],delta_old,do09_deltavir[i],do09_c200[i]),accuracy) for i in range(len(do09_m200))]
do09_mvir_p = [round(do09_mvir[i]*(do09_m200_p[i]/do09_m200[i]),accuracy) for i in range(len(do09_m200))]
do09_mvir_m = [round(do09_mvir[i]*(do09_m200_m[i]/do09_m200[i]),accuracy) for i in range(len(do09_m200))]
do09_cvir = [round(mc.Cconvert(do09_m200[i],delta_old,do09_deltavir[i],do09_c200[i]),accuracy) for i in range(len(do09_m200))]
do09_cvir_p = [round(do09_cvir[i]*(do09_c200_p[i]/do09_c200[i]),accuracy) for i in range(len(do09_m200))]
do09_cvir_m = [round(do09_cvir[i]*(do09_c200_m[i]/do09_c200[i]),accuracy) for i in range(len(do09_m200))]

for i in range(len(do09_m200)):
    print "{}, {}, {}, {}, +{}, {}, {}, +{}, {}, {}, +{}, {}, {}, +{}, {}".format(i+1,do09_clusters[i],do09_z[i],do09_c200[i],do09_c200_p[i],do09_c200_m[i],do09_m200[i],do09_m200_p[i],do09_m200_m[i],do09_cvir[i],do09_cvir_p[i],do09_cvir_m[i],do09_mvir[i],do09_mvir_p[i],do09_mvir_m[i])


