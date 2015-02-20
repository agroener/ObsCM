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
ga09_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'GA09.1']
ga09_m200 = [m200[i] for i in range(len(m200)) if short_refs[i] == 'GA09.1']
ga09_m200_p = [m200_plus[i] for i in range(len(m200_plus)) if short_refs[i] == 'GA09.1']
ga09_m200_m = [m200_minus[i] for i in range(len(m200_minus)) if short_refs[i] == 'GA09.1']
ga09_c200 = [c200[i] for i in range(len(c200)) if short_refs[i] == 'GA09.1']
ga09_c200_p = [c200_plus[i] for i in range(len(c200_plus)) if short_refs[i] == 'GA09.1']
ga09_c200_m = [c200_minus[i] for i in range(len(c200_minus)) if short_refs[i] == 'GA09.1']
ga09_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'GA09.1']
ga09_deltavir = [mc.DeltaFinder(0.3,0.7,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'GA09.1']

# testing to make sure that data are in the expected format
# make sure both concentration and mass measurements are there in all cases
assert len(ga09_m200) == len(ga09_c200)
# make sure both mass uncertainties are there in all cases, and are the same
assert len(ga09_m200_p) == len(ga09_m200), "M200 uncertainties appear to be missing for some clusters..."
assert len(ga09_m200_p) == len(ga09_m200_m), "M200 uncertainties do not have the same length..."
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(ga09_c200_p) == len(ga09_c200), "c200 uncertainties appear to be missing for some clusters..."
assert len(ga09_c200_p) == len(ga09_c200_m), "c200 uncertainties do not have the same length..."
# checking redshifts and deltavirs
assert len(ga09_z) == len(ga09_m200)
assert len(ga09_deltavir) == len(ga09_m200)

# time to convert
# propagate uncertainties by conserving fractional uncertainties
delta_orig = 200 # constant for all clusters
accuracy = 1 # number of places after the decimal to round to
ga09_mvir = [round(mc.Mconvert(ga09_m200[i],delta_orig,ga09_deltavir[i],ga09_c200[i]),accuracy) for i in range(len(ga09_m200))]
ga09_mvir_p = [round(ga09_mvir[i]*(ga09_m200_p[i]/ga09_m200[i]),accuracy) for i in range(len(ga09_m200))]
ga09_mvir_m = [round(ga09_mvir[i]*(ga09_m200_m[i]/ga09_m200[i]),accuracy) for i in range(len(ga09_m200))]
ga09_cvir = [round(mc.Cconvert(ga09_m200[i],delta_orig,ga09_deltavir[i],ga09_c200[i]),accuracy) for i in range(len(ga09_m200))]
ga09_cvir_p = [round(ga09_cvir[i]*(ga09_c200_p[i]/ga09_c200[i]),accuracy) for i in range(len(ga09_m200))]
ga09_cvir_m = [round(ga09_cvir[i]*(ga09_c200_m[i]/ga09_c200[i]),accuracy) for i in range(len(ga09_m200))]


for i in range(len(ga09_mvir)):
    print "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(i+1,ga09_clusters[i],ga09_z[i],ga09_c200[i],ga09_c200_p[i],ga09_c200_m[i],ga09_m200[i],ga09_m200_p[i],ga09_m200_m[i],ga09_cvir[i],ga09_cvir_p[i],ga09_cvir_m[i],ga09_mvir[i],ga09_mvir_p[i],ga09_mvir_m[i])
