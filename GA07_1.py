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
# or asymmetric uncertainties here to deal with
ga07_clusters = [clusters[i] for i in range(len(clusters)) if short_refs[i] == 'GA07.1']
ga07_mvir = [mvir[i] for i in range(len(mvir)) if short_refs[i] == 'GA07.1']
ga07_mvir_p = [mvir_plus[i] for i in range(len(mvir_plus)) if short_refs[i] == 'GA07.1']
ga07_mvir_m = [mvir_minus[i] for i in range(len(mvir_minus)) if short_refs[i] == 'GA07.1']
ga07_cvir = [cvir[i] for i in range(len(cvir)) if short_refs[i] == 'GA07.1']
ga07_cvir_p = [cvir_plus[i] for i in range(len(cvir_plus)) if short_refs[i] == 'GA07.1']
ga07_cvir_m = [cvir_minus[i] for i in range(len(cvir_minus)) if short_refs[i] == 'GA07.1']
ga07_z = [redshift[i] for i in range(len(redshift)) if short_refs[i] == 'GA07.1']
ga07_deltavir = [mc.DeltaFinder(0.3,0.7,redshift[i]) for i in range(len(redshift)) if short_refs[i] == 'GA07.1']

# testing to make sure that data are in the expected format
# make sure both concentration and mass measurements are there in all cases
assert len(ga07_mvir) == len(ga07_cvir)
# make sure both mass uncertainties are there in all cases, and are the same
assert len(ga07_mvir_p) == len(ga07_mvir), "Mvir uncertainties appear to be missing for some clusters..."
assert len(ga07_mvir_p) == len(ga07_mvir_m), "Mvir uncertainties do not have the same length..."
for i in range(len(ga07_mvir)):
    assert abs(ga07_mvir_p[i]) == abs(ga07_mvir_m[i]), "Mvir uncertainties are different for: {}".format(clusters[i])
# make sure both concentration uncertainties are there in all cases, and are the same
assert len(ga07_cvir_p) == len(ga07_cvir), "cvir uncertainties appear to be missing for some clusters..."
assert len(ga07_cvir_p) == len(ga07_cvir_m), "cvir uncertainties do not have the same length..."
for i in range(len(ga07_cvir)):
    assert abs(ga07_cvir_p[i]) == abs(ga07_cvir_m[i]), "cvir uncertainties are different for: {}".format(clusters[i])

# time to convert
delta_new = 200 # constant for all clusters
c_accuracy = 1 # number of places after the decimal to round to for concentration
m_accuracy = 2 # number of places after the decimal to round to for mass
ga07_m200 = [round(mc.Mconvert(ga07_mvir[i],ga07_deltavir[i],delta_new,ga07_cvir[i]),m_accuracy) for i in range(len(ga07_mvir))]
ga07_m200_err = [round(un.propagate_mass_uncertainty_notindependent(ga07_mvir[i],ga07_mvir_p[i],ga07_cvir[i],ga07_cvir_p[i],ga07_deltavir[i],delta_new),m_accuracy)
                 for i in range(len(ga07_mvir))]
ga07_c200 = [round(mc.Cconvert(ga07_mvir[i],ga07_deltavir[i],delta_new,ga07_cvir[i]),c_accuracy) for i in range(len(ga07_mvir))]
ga07_c200_err = [round(un.propagate_conc_uncertainty(ga07_mvir[i],ga07_mvir_p[i],ga07_cvir[i],ga07_cvir_p[i],ga07_deltavir[i],delta_new),c_accuracy)
                 for i in range(len(ga07_mvir))]

for i in range(len(ga07_mvir)):
    print "{}, {}, {}, {}, +/-{}, {}, +/-{}, {}, +/-{}, {}, +/-{}".format(i+1,ga07_clusters[i],ga07_z[i],ga07_c200[i],ga07_c200_err[i],ga07_m200[i],ga07_m200_err[i],ga07_cvir[i],ga07_cvir_p[i],ga07_mvir[i],ga07_mvir_p[i])
