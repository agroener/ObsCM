from xlrd import open_workbook

from astropy.table import Table
from astropy.io import ascii

import math
import matplotlib.pyplot as plt
import ipdb

import DataHandler as DH

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
z_out,methods_out = ([],[])
for i in clusters_gte10:
    tmp_mvir,tmp_mvir_p,tmp_mvir_m,tmp_cvir,tmp_cvir_p,tmp_cvir_m = ([],[],[],[],[],[])
    tmp_z,tmp_methods = ([],[])
    for j in range(len(clusters)):
        if clusters[j] == i:
            if math.isnan(float(mvir[j])) or math.isnan(float(cvir[j])):
                continue
            else:
                tmp_mvir.append(mvir[j])
                tmp_mvir_p.append(mvir_plus[j])
                tmp_mvir_m.append(mvir_minus[j])
                tmp_cvir.append(cvir[j])
                tmp_cvir_p.append(cvir_plus[j])
                tmp_cvir_m.append(cvir_minus[j])
                tmp_z.append(redshift[j])
                tmp_methods.append(methods[j])
    mvir_out.append(tmp_mvir)
    mvir_p_out.append(tmp_mvir_p)
    mvir_m_out.append(tmp_mvir_m)
    cvir_out.append(tmp_cvir)
    cvir_p_out.append(tmp_cvir_p)
    cvir_m_out.append(tmp_cvir_m)
    z_out.append(tmp_z)
    methods_out.append(tmp_methods)
    
nplots = len(clusters_gte10) # subplots do not auto-update; need to do manually
f, axes = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
f.set_size_inches(7,7,forward=True)
# loop structure only works for 3 clusters in a 2x2 grid format
for i,j in enumerate(clusters_gte10):
    if i <= 1:
        #axes[0][i].set_title('{}'.format(j))
        axes[0][i].set_xlim(0,175)
        axes[0][i].set_ylim(0,30)
    else:
        #axes[1][0].set_title('{}'.format(j))
        axes[1][0].set_xlim(0,175)
        axes[1][0].set_ylim(0,30)
    print j
    for k in range(len(mvir_out[i])):
        print "Measurement number {}".format(k+1)
        flag = DH.scatter_flag(k,mvir_out[i],mvir_p_out[i],mvir_m_out[i],
                               cvir_out[i],cvir_p_out[i],cvir_m_out[i])
        col = DH.colorselect(methods_out[i][k])
        if i <= 1:
            DH.plotwithflag(axes[0][i],flag,col,mvir_out[i][k],mvir_p_out[i][k],mvir_m_out[i][k],
                        cvir_out[i][k],cvir_p_out[i][k],cvir_m_out[i][k],z_out[i][k],witherrors=True)
        else:
            DH.plotwithflag(axes[1][0],flag,col,mvir_out[i][k],mvir_p_out[i][k],mvir_m_out[i][k],
                        cvir_out[i][k],cvir_p_out[i][k],cvir_m_out[i][k],z_out[i][k],witherrors=True)

axes[1][1].scatter(1e6,1e6,color=DH.colorselect('X-ray'),label='X-ray')
axes[1][1].scatter(1e6,1e6,color=DH.colorselect('WL'),label='WL')
axes[1][1].scatter(1e6,1e6,color=DH.colorselect('SL'),label='SL')
axes[1][1].scatter(1e6,1e6,color=DH.colorselect('WL+SL'),label='WL+SL')
axes[1][1].scatter(1e6,1e6,color=DH.colorselect('CM'),label='CM')
axes[1][1].scatter(1e6,1e6,color=DH.colorselect('LOSVD'),label='LOSVD')
axes[1][1].set_xlim(0,175)
axes[1][1].set_ylim(0,30)

# over-plot cm relations
# first CO07_1
mlistCO07_1,clistCO07_1,zCO07_1 = DH.cmrelation_co07_1(1e13,1e16,0)
mlistCO07_1 = [mlistCO07_1[i]/1.e13 for i in range(len(mlistCO07_1))]
mlistCO07_1_p,clistCO07_1_p,zCO07_1 = DH.cmrelation_co07_1(1e13,1e16,0,c0=20.9,alpha=-0.02)
mlistCO07_1_p = [mlistCO07_1_p[i]/1.e13 for i in range(len(mlistCO07_1_p))]
mlistCO07_1_m,clistCO07_1_m,zCO07_1 = DH.cmrelation_co07_1(1e13,1e16,0,c0=8.7,alpha=-0.26)
mlistCO07_1_m = [mlistCO07_1_m[i]/1.e13 for i in range(len(mlistCO07_1_m))]

mlistBU01_1,clistBU01_1,zBU01_1 = DH.cmrelation_bu01_1(1e13,1e16,0)
mlistBU01_1 = [mlistBU01_1[i]/1.e13 for i in range(len(mlistBU01_1))]

mlistHE07_1,clistHE07_1,zHE07_1 = DH.cmrelation_he07_1(1e13,1e16,0)
mlistHE07_1 = [mlistHE07_1[i]/1.e13 for i in range(len(mlistHE07_1))]

axes[0][0].plot(mlistCO07_1,clistCO07_1,color='b',linewidth=2,linestyle='--')
axes[0][0].plot(mlistBU01_1,clistBU01_1,color='y',linewidth=2,linestyle='--')
axes[0][0].plot(mlistHE07_1,clistHE07_1,color='orange',linewidth=2,linestyle='--')
axes[0][0].fill_between(mlistCO07_1,clistCO07_1_p,clistCO07_1_m,alpha=0.25)

axes[0][1].plot(mlistCO07_1,clistCO07_1,color='b',linewidth=2,linestyle='--')
axes[0][1].plot(mlistBU01_1,clistBU01_1,color='y',linewidth=2,linestyle='--')
axes[0][1].plot(mlistHE07_1,clistHE07_1,color='orange',linewidth=2,linestyle='--')
axes[0][1].fill_between(mlistCO07_1,clistCO07_1_p,clistCO07_1_m,alpha=0.25)

axes[1][0].plot(mlistCO07_1,clistCO07_1,color='b',linewidth=2,linestyle='--')
axes[1][0].plot(mlistBU01_1,clistBU01_1,color='y',linewidth=2,linestyle='--')
axes[1][0].plot(mlistHE07_1,clistHE07_1,color='orange',linewidth=2,linestyle='--')
axes[1][0].fill_between(mlistCO07_1,clistCO07_1_p,clistCO07_1_m,alpha=0.25)

f.subplots_adjust(hspace=0)
f.subplots_adjust(wspace=0)

plt.legend(loc=0,scatterpoints=1)
plt.show()
