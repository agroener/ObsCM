from xlrd import open_workbook

from astropy.table import Table
from astropy.io import ascii

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

# clusters with 10 or more measurements
clusters_gte10 = set([i for i in clusters if clusters.count(i) >= 10])
mvir_out,mvir_p_out,mvir_m_out,cvir_out,cvir_p_out,cvir_m_out = ([],[],[],[],[],[])
for i in clusters_gte10:
    tmp_mvir,tmp_mvir_p,tmp_mvir_m,tmp_cvir,tmp_cvir_p,tmp_cvir_m = ([],[],[],[],[],[])
    for j in range(len(clusters)):
        if clusters[j] == i:
            tmp_mvir.append(mvir[j])
            tmp_mvir_p.append(mvir_plus[j])
            tmp_mvir_m.append(mvir_minus[j])
            tmp_cvir.append(cvir[j])
            tmp_cvir_p.append(cvir_plus[j])
            tmp_cvir_m.append(cvir_minus[j])
    mvir_out.append(tmp_mvir)
    mvir_p_out.append(tmp_mvir_p)
    mvir_m_out.append(tmp_mvir_m)
    cvir_out.append(tmp_cvir)
    cvir_p_out.append(tmp_cvir_p)
    cvir_m_out.append(tmp_cvir_m)
ipdb.set_trace()

