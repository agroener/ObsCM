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
og13_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'OG13.1']
og13_mvir = [mvir[i] for i in range(len(mvir)) if short_refs[i] == 'OG13.1']
og13_mvir_p = [mvir_plus[i] for i in range(len(mvir_plus)) if short_refs[i] == 'OG13.1']
og13_mvir_m = [mvir_minus[i] for i in range(len(mvir_minus)) if short_refs[i] == 'OG13.1']
og13_cvir = [cvir[i] for i in range(len(cvir)) if short_refs[i] == 'OG13.1']
og13_cvir_p = [cvir_plus[i] for i in range(len(cvir_plus)) if short_refs[i] == 'OG13.1']
og13_cvir_m = [cvir_minus[i] for i in range(len(cvir_minus)) if short_refs[i] == 'OG13.1']
og13_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'OG13.1']
og13_deltavir = [mc.DeltaFinder(0.275,0.725,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'OG13.1']

# testing to make sure that data are in the expected format
# a few orders of business first
assert len(og13_clusters) == len(og13_z)
assert len(og13_z) == len(og13_deltavir)
# make sure both concentration and mass measurements are there in all cases
assert len(og13_mvir) == len(og13_cvir)
assert 'TBD' not in og13_mvir
assert 'TBD' not in og13_cvir
# make sure both mass uncertainties are there in all cases, and are the same
assert len(og13_mvir_p) == len(og13_mvir), "Mvir uncertainties appear to be missing for some clusters..."
assert len(og13_mvir_p) == len(og13_mvir_m), "Mvir uncertainties do not have the same length..."
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(og13_cvir_p) == len(og13_cvir), "cvir uncertainties appear to be missing for some clusters..."
assert len(og13_cvir_p) == len(og13_cvir_m), "cvir uncertainties do not have the same length..."

# time to convert; for asymmetric uncertainties, I've decided to keep the fractional error constant
delta_new = 200 # constant for all clusters
accuracy = 2 # number of places after the decimal to round to
og13_m200 = [round(mc.Mconvert(og13_mvir[i],og13_deltavir[i],delta_new,og13_cvir[i]),accuracy) for i in range(len(og13_mvir))]
og13_m200_p = [round(og13_m200[i]*(og13_mvir_p[i]/og13_mvir[i]),accuracy) for i in range(len(og13_mvir))]
og13_m200_m = [round(og13_m200[i]*(og13_mvir_m[i]/og13_mvir[i]),accuracy) for i in range(len(og13_mvir))]
og13_c200 = [round(mc.Cconvert(og13_mvir[i],og13_deltavir[i],delta_new,og13_cvir[i]),accuracy) for i in range(len(og13_mvir))]
og13_c200_p = [round(og13_c200[i]*(og13_cvir_p[i]/og13_cvir[i]),accuracy) for i in range(len(og13_mvir))]
og13_c200_m = [round(og13_c200[i]*(og13_cvir_m[i]/og13_cvir[i]),accuracy) for i in range(len(og13_mvir))]

for i in range(len(og13_mvir)):
    print "{}, {}, {}, {}, +{}, {}, {}, +{}, {}, {}, +{}, {}, {}, +{}, {}".format(i+1,og13_clusters[i],og13_z[i],og13_c200[i],og13_c200_p[i],og13_c200_m[i],og13_m200[i],og13_m200_p[i],og13_m200_m[i],og13_cvir[i],og13_cvir_p[i],og13_cvir_m[i],og13_mvir[i],og13_mvir_p[i],og13_mvir_m[i])



