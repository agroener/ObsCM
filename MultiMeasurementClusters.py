from xlrd import open_workbook

from astropy.table import Table
from astropy.io import ascii

import math
from matplotlib import pyplot
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

def clustersgte10():
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
            axes[0][i].set_xlim(1e14,2e16)
            axes[0][i].set_ylim(1,30)
            axes[0][i].text(3e14,1.5,'{}'.format(j))
            axes[0][i].text(3e14,1.25,'z={}'.format(z_out[i][0]))
            axes[0][i].set_xscale('log')
            axes[0][i].set_yscale('log')
        else:
            axes[1][0].set_xlim(1e14,2e16)
            axes[1][0].set_ylim(1,30)
            axes[1][0].text(3e14,1.5,'{}'.format(j))
            axes[1][0].text(3e14,1.25,'z={}'.format(z_out[i][0]))
            axes[1][0].set_xscale('log')
            axes[1][0].set_yscale('log')
        print j
        for k in range(len(mvir_out[i])):
            print "Measurement number {}".format(k+1)
            flag = DH.scatter_flag(k,mvir_out[i],mvir_p_out[i],mvir_m_out[i],
                                   cvir_out[i],cvir_p_out[i],cvir_m_out[i])
            col = DH.colorselect(methods_out[i][k])
            if i <= 1:
                DH.plotwithflag(axes[0][i],flag,col,mvir_out[i][k]*1e14,mvir_p_out[i][k]*1e14,mvir_m_out[i][k]*1e14,
                                cvir_out[i][k],cvir_p_out[i][k],cvir_m_out[i][k],z_out[i][k],witherrors=True)
            else:
                DH.plotwithflag(axes[1][0],flag,col,mvir_out[i][k]*1e14,mvir_p_out[i][k]*1e14,mvir_m_out[i][k]*1e14,
                                cvir_out[i][k],cvir_p_out[i][k],cvir_m_out[i][k],z_out[i][k],witherrors=True)

    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('X-ray'),label='X-ray')
    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('WL'),label='WL')
    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('SL'),label='SL')
    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('WL+SL'),label='WL+SL')
    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('CM'),label='CM')
    axes[1][1].scatter(1e6,1e6,color=DH.colorselect('LOSVD'),label='LOSVD')
    axes[1][1].set_xlim(1e14,2e16)
    axes[1][1].set_ylim(1,30)

    # over-plot cm relations
    # CO07_1
    mlistCO07_1_1689,clistCO07_1_1689,zCO07_1_1689 = DH.cmrelation_co07_1(1e14,2e16,z_out[0][0])
    mlistCO07_1_p_1689,clistCO07_1_p_1689,zCO07_1_1689 = DH.cmrelation_co07_1(1e14,2e16,z_out[0][0],c0=20.9,alpha=-0.02)
    mlistCO07_1_m_1689,clistCO07_1_m_1689,zCO07_1_1689 = DH.cmrelation_co07_1(1e14,2e16,z_out[0][0],c0=8.7,alpha=-0.26)
    mlistCO07_1_2137,clistCO07_1_2137,zCO07_1_2137 = DH.cmrelation_co07_1(1e14,2e16,z_out[1][0])
    mlistCO07_1_p_2137,clistCO07_1_p_2137,zCO07_1_2137 = DH.cmrelation_co07_1(1e14,2e16,z_out[1][0],c0=20.9,alpha=-0.02)
    mlistCO07_1_m_2137,clistCO07_1_m_2137,zCO07_1_2137 = DH.cmrelation_co07_1(1e14,2e16,z_out[1][0],c0=8.7,alpha=-0.26)
    mlistCO07_1_1835,clistCO07_1_1835,zCO07_1_1835 = DH.cmrelation_co07_1(1e14,2e16,z_out[2][0])
    mlistCO07_1_p_1835,clistCO07_1_p_1835,zCO07_1_1835 = DH.cmrelation_co07_1(1e14,2e16,z_out[2][0],c0=20.9,alpha=-0.02)
    mlistCO07_1_m_1835,clistCO07_1_m_1835,zCO07_1_1835 = DH.cmrelation_co07_1(1e14,2e16,z_out[2][0],c0=8.7,alpha=-0.26)
    # GR14_1
    mlistGR14_1_1689,clistGR14_1_1689,zGR14_1_1689 = DH.cmrelation_gr14_1(1e14,2e16,z_out[0][0])
    mlistGR14_1_p_1689,clistGR14_1_p_1689,zGR14_1_1689 = DH.cmrelation_gr14_1(1e14,2e16,z_out[0][0],c0=4.797,alpha=-0.049)
    mlistGR14_1_m_1689,clistGR14_1_m_1689,zGR14_1_1689 = DH.cmrelation_gr14_1(1e14,2e16,z_out[0][0],c0=4.753,alpha=-0.063)
    mlistGR14_1_2137,clistGR14_1_2137,zGR14_1_2137 = DH.cmrelation_gr14_1(1e14,2e16,z_out[1][0])
    mlistGR14_1_p_2137,clistGR14_1_p_2137,zGR14_1_2137 = DH.cmrelation_gr14_1(1e14,2e16,z_out[1][0],c0=4.797,alpha=-0.049)
    mlistGR14_1_m_2137,clistGR14_1_m_2137,zGR14_1_2137 = DH.cmrelation_gr14_1(1e14,2e16,z_out[1][0],c0=4.753,alpha=-0.063)
    mlistGR14_1_1835,clistGR14_1_1835,zGR14_1_1835 = DH.cmrelation_gr14_1(1e14,2e16,z_out[2][0])
    mlistGR14_1_p_1835,clistGR14_1_p_1835,zGR14_1_1835 = DH.cmrelation_gr14_1(1e14,2e16,z_out[2][0],c0=4.797,alpha=-0.049)
    mlistGR14_1_m_1835,clistGR14_1_m_1835,zGR14_1_1835 = DH.cmrelation_gr14_1(1e14,2e16,z_out[2][0],c0=4.753,alpha=-0.063)
    # BU01_1
    mlistBU01_1_1689,clistBU01_1_1689,zBU01_1_1689 = DH.cmrelation_bu01_1(1e14,2e16,z_out[0][0])
    mlistBU01_1_2137,clistBU01_1_2137,zBU01_1_2137 = DH.cmrelation_bu01_1(1e14,2e16,z_out[1][0])
    mlistBU01_1_1835,clistBU01_1_1835,zBU01_1_1835 = DH.cmrelation_bu01_1(1e14,2e16,z_out[2][0])
    # HE07_1
    mlistHE07_1_1689,clistHE07_1_1689,zHE07_1_1689 = DH.cmrelation_he07_1(1e14,2e16,z_out[0][0])
    mlistHE07_1_2137,clistHE07_1_2137,zHE07_1_2137 = DH.cmrelation_he07_1(1e14,2e16,z_out[1][0])
    mlistHE07_1_1835,clistHE07_1_1835,zHE07_1_1835 = DH.cmrelation_he07_1(1e14,2e16,z_out[2][0])
    # PR11_1
    #mlistPR11_1_1689,clistPR11_1_1689,zPR11_1_1689 = DH.cmrelation_pr11_1(1e14,2e16,z_out[0][0])
    #mlistPR11_1_2137,clistPR11_1_2137,zPR11_1_2137 = DH.cmrelation_pr11_1(1e14,2e16,z_out[1][0])
    #mlistPR11_1_1835,clistPR11_1_1835,zPR11_1_1835 = DH.cmrelation_pr11_1(1e14,2e16,z_out[2][0])

    # Abell 1689
    axes[0][0].plot(mlistCO07_1_1689,clistCO07_1_1689,color='b',linewidth=2,linestyle='--')
    #axes[0][0].plot(mlistBU01_1_1689,clistBU01_1_1689,color='y',linewidth=2,linestyle='--')
    #axes[0][0].plot(mlistHE07_1_1689,clistHE07_1_1689,color='orange',linewidth=2,linestyle='--')
    axes[0][0].plot(mlistGR14_1_1689,clistGR14_1_1689,color='green',linewidth=2,linestyle='--')
    #axes[0][0].plot(mlistPR11_1_1689,clistPR11_1_1689,color='green',linewidth=2,linestyle='--')
    axes[0][0].fill_between(mlistCO07_1_1689,clistCO07_1_p_1689,clistCO07_1_m_1689,alpha=0.25,color='blue')
    axes[0][0].fill_between(mlistGR14_1_1689,clistGR14_1_p_1689,clistGR14_1_m_1689,alpha=0.25,color='green')
    # MS 2137
    axes[0][1].plot(mlistCO07_1_2137,clistCO07_1_2137,color='b',linewidth=2,linestyle='--')
    #axes[0][1].plot(mlistBU01_1_2137,clistBU01_1_2137,color='y',linewidth=2,linestyle='--')
    #axes[0][1].plot(mlistHE07_1_2137,clistHE07_1_2137,color='orange',linewidth=2,linestyle='--')
    axes[0][1].plot(mlistGR14_1_2137,clistGR14_1_2137,color='green',linewidth=2,linestyle='--')
    #axes[0][1].plot(mlistPR11_1_2137,clistPR11_1_2137,color='green',linewidth=2,linestyle='--')
    axes[0][1].fill_between(mlistCO07_1_2137,clistCO07_1_p_2137,clistCO07_1_m_2137,alpha=0.25,color='blue')
    axes[0][1].fill_between(mlistGR14_1_2137,clistGR14_1_p_2137,clistGR14_1_m_2137,alpha=0.25,color='green')
    # Abell 1835
    axes[1][0].plot(mlistCO07_1_1835,clistCO07_1_1835,color='b',linewidth=2,linestyle='--')
    #axes[1][0].plot(mlistBU01_1_1835,clistBU01_1_1835,color='y',linewidth=2,linestyle='--')
    #axes[1][0].plot(mlistHE07_1_1835,clistHE07_1_1835,color='orange',linewidth=2,linestyle='--')
    axes[1][0].plot(mlistGR14_1_1835,clistGR14_1_1835,color='green',linewidth=2,linestyle='--')
    #axes[1][0].plot(mlistPR11_1_1835,clistPR11_1_1835,color='green',linewidth=2,linestyle='--')
    axes[1][0].fill_between(mlistCO07_1_1835,clistCO07_1_p_1835,clistCO07_1_m_1835,alpha=0.25,color='blue')
    axes[1][0].fill_between(mlistGR14_1_1835,clistGR14_1_p_1835,clistGR14_1_m_1835,alpha=0.25,color='green')

    axes[1][1].plot(1e6,1e6,color='b',linewidth=2,linestyle='--',label='CO07')
    #axes[1][1].plot(1e6,1e6,color='y',linewidth=2,linestyle='--',label='BU01')
    #axes[1][1].plot(1e6,1e6,color='orange',linewidth=2,linestyle='--',label='HE07')
    axes[1][1].plot(1e6,1e6,color='green',linewidth=2,linestyle='--',label='GR14')
    #axes[1][1].plot(1e6,1e6,color='green',linewidth=2,linestyle='--',label='PR11')
    axes[1][1].set_xscale('log')
    axes[1][1].set_yscale('log')

    f.subplots_adjust(hspace=0)
    f.subplots_adjust(wspace=0)
    f.text(0.5, 0.04, r'$\mathrm{M_{vir} (M_{\odot})}$', ha='center', va='center', fontsize=18)
    f.text(0.06, 0.5, r'$\mathrm{c_{vir} (1+z)}$', ha='center', va='center', rotation='vertical', fontsize=18)

    plt.legend(loc=0,scatterpoints=1,frameon=False)
    plt.show()


if __name__ == '__main__':
    
    clustersgte10()
