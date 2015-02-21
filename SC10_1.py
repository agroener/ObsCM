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
sc10_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'SC10.1']
sc10_m200 = [m200[i] for i in range(len(m200)) if short_refs[i] == 'SC10.1']
sc10_m200_p = [m200_plus[i] for i in range(len(m200_plus)) if short_refs[i] == 'SC10.1']
sc10_m200_m = [m200_minus[i] for i in range(len(m200_minus)) if short_refs[i] == 'SC10.1']
sc10_c200 = [c200[i] for i in range(len(c200)) if short_refs[i] == 'SC10.1']
sc10_c200_p = [c200_plus[i] for i in range(len(c200_plus)) if short_refs[i] == 'SC10.1']
sc10_c200_m = [c200_minus[i] for i in range(len(c200_minus)) if short_refs[i] == 'SC10.1']
sc10_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'SC10.1']
sc10_deltavir = [mc.DeltaFinder(0.27,0.73,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'SC10.1']

# testing to make sure that data are in the expected format
# make sure both concentration and mass measurements are there in all cases
assert len(sc10_m200) == len(sc10_c200)
# make sure both mass uncertainties are there in all cases, and are the same
assert len(sc10_m200_p) == len(sc10_m200), "M200 uncertainties appear to be missing for some clusters..."
assert len(sc10_m200_p) == len(sc10_m200_m), "M200 uncertainties do not have the same length..."
for i in range(len(sc10_m200)):
    assert abs(sc10_m200_p[i]) == abs(sc10_m200_m[i]), "Mass uncertainties differ for {}".format(sc10_clusters[i])
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(sc10_c200_p) == len(sc10_c200), "c200 uncertainties appear to be missing for some clusters..."
assert len(sc10_c200_p) == len(sc10_c200_m), "c200 uncertainties do not have the same length..."
for i in range(len(sc10_m200)):
    assert abs(sc10_c200_p[i]) == abs(sc10_c200_m[i]), "Conc. uncertainties differ for {}".format(sc10_clusters[i])
# checking redshifts and deltavirs
assert len(sc10_z) == len(sc10_m200)
assert len(sc10_deltavir) == len(sc10_m200)

# time to convert
# propagate uncertainties by conserving fractional uncertainties
# decided to do this for all CO09 clusters because one or more had asymmetric uncertainties
delta_orig = 200 # constant for all clusters
accuracy = 2 # number of places after the decimal to round to
sc10_mvir = [round(mc.Mconvert(sc10_m200[i],delta_orig,sc10_deltavir[i],sc10_c200[i]),accuracy) for i in range(len(sc10_m200))]
sc10_mvir_err = [round(un.propagate_mass_uncertainty_notindependent(sc10_m200[i],sc10_m200_p[i],sc10_c200[i],sc10_c200_p[i],delta_orig,sc10_deltavir[i]),accuracy) for i in range(len(sc10_m200))]
sc10_cvir = [round(mc.Cconvert(sc10_m200[i],delta_orig,sc10_deltavir[i],sc10_c200[i]),accuracy) for i in range(len(sc10_m200))]
sc10_cvir_err = [round(un.propagate_conc_uncertainty(sc10_m200[i],sc10_m200_p[i],sc10_c200[i],sc10_c200_p[i],delta_orig,sc10_deltavir[i]),accuracy) for i in range(len(sc10_m200))]

for i in range(len(sc10_mvir)):
    print "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, -{}, {}, {}, -{}".format(i+1,sc10_clusters[i],sc10_z[i],sc10_c200[i],sc10_c200_p[i],sc10_c200_m[i],sc10_m200[i],sc10_m200_p[i],sc10_m200_m[i],sc10_cvir[i],sc10_cvir_err[i],sc10_cvir_err[i],sc10_mvir[i],sc10_mvir_err[i],sc10_mvir_err[i])

