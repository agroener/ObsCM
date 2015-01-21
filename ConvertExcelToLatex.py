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

# Create astropy.Table object

## Other things to check: 'infty' turns into infinity symbol
## Need to compress errors as superscript and subscript of measurements of c/M to save space.

t = Table()

clusters = [i.encode('utf-8') for i in clusters]
t['Clusters'] = clusters
t['Redshift'] = redshift
t['Method'] = methods
t['c200'] = c200
t['c200+err'] = c200_plus
t['c200-err'] = c200_minus
t['M200'] = m200
t['M200+err'] = m200_plus
t['M200-err'] = m200_minus
t['cvir'] = cvir
t['cvir+err'] = cvir_plus
t['cvir-err'] = cvir_minus
t['Mvir'] = mvir
t['Mvir+err'] = mvir_plus
t['Mvir-err'] = mvir_minus
t['Ref.'] = short_refs
t['Orig. Convention'] = orig_convention
t['Cosmology'] = cosmology

    
# Write data out to latex file
t.write('cm_data_test.tex', format='latex')
