import numpy as np
from xlrd import open_workbook
import matplotlib.pyplot as plt

import ipdb

## The general conversion process will take data in one convention (either 200 or virial), and
## convert these measurements to a common cosmology and uncertainty convention.

## Functions for converting between overdensity convention (if needed) come from MConvert_Personal.py.
import MConvert_Personal as MC
## Functions for converting between cosmologies come from CosmoConvert.py.
import CosmoConvert as CC
## Functions for converting confidence intervals to symmetric uncertainties come from Uncertainties.py.
import Uncertainties as UN

# Get the data form the excel sheet cm_data.xlsx
import DataHandler as DH
clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology=DH.startup()

# Take all virial measurements and normalize uncertainties
mvir_norm,mvir_p_norm,mvir_m_norm,cvir_norm,cvir_p_norm,cvir_m_norm = ([],[],[],[],[],[])
methods_norm,z_norm = ([],[])
for i in range(len(clusters)):
    if mvir[i] not in [u'TBD',u'nan']:
        tmp = UN.normalize_uncertainty(mvir[i],mvir_plus[i],mvir_minus[i],cvir[i],cvir_plus[i],cvir_minus[i])
        mvir_norm.append(tmp[0])
        mvir_p_norm.append(tmp[1])
        mvir_m_norm.append(tmp[2])
        cvir_norm.append(tmp[3])
        cvir_p_norm.append(tmp[4])
        cvir_m_norm.append(tmp[5])
        methods_norm.append(methods[i])
        z_norm.append(redshift[i])
        
# Plotting the results
plt.figure(figsize=(8,8))
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\mathrm{M_{vir}/M_{\odot}}$')
plt.ylabel(r'$\mathrm{c_{vir}}$')
for i in range(len(mvir_norm)):
    if methods_norm[i] in ['WL']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='purple')
    if methods_norm[i] in ['SL']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='red')
    if methods_norm[i] in ['WL+SL']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='black')
    if methods_norm[i] in ['X-ray']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='green')
    if methods_norm[i] in ['CM']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='blue')
    if methods_norm[i] in ['LOSVD']:
        plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]],xerr=[mvir_p_norm[i]*1e14],color='orange')
plt.show()

