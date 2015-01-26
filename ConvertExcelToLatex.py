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
for i in range(len(c200)):
    if c200[i] == u'TBD':
        c200[i] = '\mathrm{TBD}'
    if c200_plus[i] == u'TBD':
        c200_plus[i] = '\mathrm{TBD}'
    if c200_minus[i] == u'TBD':
        c200_minus[i] = '\mathrm{TBD}'
    if m200[i] == u'TBD':
        m200[i] = '\mathrm{TBD}'
    if m200_plus[i] == u'TBD':
        m200_plus[i] = '\mathrm{TBD}'
    if m200_minus[i] == u'TBD':
        m200_minus[i] = '\mathrm{TBD}'
    if cvir[i] == u'TBD':
        cvir[i] = '\mathrm{TBD}'
    if cvir_plus[i] == u'TBD':
        cvir_plus[i] = '\mathrm{TBD}'
    if cvir_minus[i] == u'TBD':
        cvir_minus[i] = '\mathrm{TBD}'
    if mvir[i] == u'TBD':
        mvir[i] = '\mathrm{TBD}'
    if mvir_plus[i] == u'TBD':
        mvir_plus[i] = '\mathrm{TBD}'
    if mvir_minus[i] == u'TBD':
        mvir_minus[i] = '\mathrm{TBD}'
for i in range(len(c200)):
    if c200[i] == u'nan':
        c200[i] = '-'
    if c200_plus[i] == u'nan':
        c200_plus[i] = ''
    if c200_minus[i] == u'nan':
        c200_minus[i] = ''
    if m200[i] == u'nan':
        m200[i] = '-'
    if m200_plus[i] == u'nan':
        m200_plus[i] = ''
    if m200_minus[i] == u'nan':
        m200_minus[i] = ''
    if cvir[i] == u'nan':
        cvir[i] = '-'
    if cvir_plus[i] == u'nan':
        cvir_plus[i] = ''
    if cvir_minus[i] == u'nan':
        cvir_minus[i] = ''
    if mvir[i] == u'nan':
        mvir[i] = '-'
    if mvir_plus[i] == u'nan':
        mvir_plus[i] = ''
    if mvir_minus[i] == u'nan':
        mvir_minus[i] = ''
for i in range(len(c200)):
    if c200_plus[i] not in ['', '\mathrm{TBD}']:
        c200_plus[i] = "+{}".format(c200_plus[i])
    if m200_plus[i] not in ['', '\mathrm{TBD}']:
        m200_plus[i] = "+{}".format(m200_plus[i])
    if cvir_plus[i] not in ['', '\mathrm{TBD}']:
        cvir_plus[i] = "+{}".format(cvir_plus[i])
    if mvir_plus[i] not in ['', '\mathrm{TBD}']:
        mvir_plus[i] = "+{}".format(mvir_plus[i])
for i in range(len(c200)):
    if c200_plus[i] == '+infty':
        c200_plus[i] = "+\infty"
    if m200_plus[i] == '+infty':
        m200_plus[i] = "+\infty"
    if cvir_plus[i] == '+infty':
        cvir_plus[i] = "+\infty"
    if mvir_plus[i] == '+infty':
        mvir_plus[i] = "+\infty"
for i in range(len(cosmology)):
    cosmology[i] = cosmology[i].strip('(').strip(')')
for i in range(len(short_refs)):
    short_refs[i] = "\citet{"+"{}".format(short_refs[i])+"}"

t = Table()

clusters = [i.encode('utf-8') for i in clusters]
t['Clusters'] = clusters
t['Redshift'] = redshift
t['Method'] = methods
c200_all = ["${"+"{}".format(c200[i])+"}^{" + "{}".format(c200_plus[i]) +"}_{" + "{}".format(c200_minus[i]) +"}$" for i in range(len(c200))]
t['c200'] = c200_all
M200_all = ["${"+"{}".format(m200[i])+"}^{" + "{}".format(m200_plus[i]) +"}_{" + "{}".format(m200_minus[i]) +"}$" for i in range(len(m200))]
t['M200'] = M200_all
cvir_all = ["${"+"{}".format(cvir[i])+"}^{" + "{}".format(cvir_plus[i]) +"}_{" + "{}".format(cvir_minus[i]) +"}$" for i in range(len(cvir))]
t['cvir'] = cvir_all
Mvir_all = ["${"+"{}".format(mvir[i])+"}^{" + "{}".format(mvir_plus[i]) +"}_{" + "{}".format(mvir_minus[i]) +"}$" for i in range(len(mvir))]
t['Mvir'] = Mvir_all
t['Ref.'] = short_refs
orig_convention_corr = []
for i in range(len(orig_convention)):
    try:
        orig_convention_corr.append(int(orig_convention[i]))
    except ValueError:
        orig_convention_corr.append(orig_convention[i])
t['Orig. Convention'] = orig_convention_corr
t['Cosmology'] = cosmology
    
# Write data out to latex file
t.write('/Users/groenera/Desktop/Dropbox/Private/Research/Papers/My_Publications/ConcentrationMassProject/Iteration1/cm_data_test.tex', format='latex')
