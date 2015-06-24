from __future__ import division
from astropy.io.ascii import read
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq
import random

import ipdb

#################
### Functions ###
#################

# Only Flat Cosmologies
def Ez(z,Omega_m_0 = 0.3, Omega_L_0 = 0.7):
    return np.sqrt(Omega_m_0*(1+z)**3 + Omega_L_0)

# Mass in units of 10^14 M_sun
def Lx(Mass, alpha = -0.13, beta = 1.60, gamma_z = -1.74, gamma_ss = 7./3., sigma = 0.11, z = 0, Omega_m_0 = 0.3, Omega_L_0 = 0.7):
    return 10**(random.normalvariate(0,sigma) + alpha + beta*np.log10(Mass)+(gamma_z+gamma_ss)*np.log10(Ez(z,Omega_m_0,Omega_L_0)))

# Functions for finding projected concentrations
def j(p,q,phi,theta):
    return np.sin(theta)**2/(p**2*q**2) + np.cos(theta)**2*np.cos(phi)**2/q**2 + np.cos(theta)**2*np.sin(phi)**2

def k(p,q,phi,theta):
    return (1/q**2 - 1/p**2)*np.sin(phi)*np.cos(phi)*np.cos(theta)

def l(p,q,phi,theta):
    return np.sin(phi)**2/q**2 + np.cos(phi)**2/p**2

def Q(p,q,phi,theta):
    return np.sqrt((j(p,q,phi,theta)+l(p,q,phi,theta)-np.sqrt((j(p,q,phi,theta)-l(p,q,phi,theta))**2+4*k(p,q,phi,theta)**2))/(j(p,q,phi,theta)+l(p,q,phi,theta)+np.sqrt((j(p,q,phi,theta)-l(p,q,phi,theta))**2+4*k(p,q,phi,theta)**2)))

def f(p,q,phi,theta):
    return np.sin(theta)**2*np.sin(phi)**2/q**2 + np.sin(theta)**2*np.cos(phi)**2/p**2 + np.cos(theta)**2

def e_delta(p,q,phi,theta):
    return np.sqrt(p*q/Q(p,q,phi,theta)) * f(p,q,phi,theta)**(3/4)

def deltaC_deltac_ratio(p,q,phi,theta):
    return 1/e_delta(p,q,phi,theta)

def delta(c):
    return (200/3) * (c**3 / (np.log(1+c) - c/(1+c)))

def residuals(C,c,q,thetaval,ProlOrObl):
    if ProlOrObl == 'prol':
        return delta(C) - deltaC_deltac_ratio(q,q,0,thetaval) * delta(c)
    elif ProlOrObl == 'obl':
        return delta(C) - deltaC_deltac_ratio(1,q,np.pi/2,thetaval) * delta(c)

def conc_finder_pro(c,thetalist,q):
    outlist = []
    for theta in thetalist:
        p0 = [c]
        outlist.append(leastsq(residuals,p0, args=(c,q,theta,'prol'))[0])
    return outlist

def startup_sims():
    try:
        #ipdb.set_trace()
        raw_data = read('/Users/groenera/Desktop/Dropbox/Private/Research/GroupMeetings/Meeting#60/concs_m200s_data.txt',delimiter=',',guess=False)
    except:
        print("Cannot find datafile...")
        return

    raw_concs = raw_data['c200']
    raw_fof = raw_data['mfof']
    raw_masses = raw_data['m200']

    high_indices = [i for i in range(len(raw_fof)) if raw_fof[i] > 1.2e14]
    h_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in high_indices]
    h_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in high_indices]

    low_indices = [i for i in range(len(raw_masses)) if raw_fof[i] < 2.6e13]
    l_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in low_indices]
    l_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in low_indices]

    med_indices = [i for i in range(len(raw_masses)) if i not in high_indices and i not in low_indices]
    m_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in med_indices]
    m_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in med_indices]

    return l_concs,l_masses,m_concs,m_masses,h_concs,h_masses

if __name__ == "__main__":
    # intrinsic 3D masses and concentrations
    l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()
    
    # projected concentrations
    l_concs_p = [conc_finder_pro(l_concs[i],[0],0.5) for i in range(len(l_concs))]
    m_concs_p = [conc_finder_pro(m_concs[i],[0],0.5) for i in range(len(m_concs))]
    h_concs_p = [conc_finder_pro(h_concs[i],[0],0.5) for i in range(len(h_concs))]

    # Luminosities based on intrinsic masses
    Lx_l = [Lx(i/1.e14)*(1.e44) for i in l_masses]
    Lx_m = [Lx(i/1.e14)*(1.e44) for i in m_masses]
    Lx_h = [Lx(i/1.e14)*(1.e44) for i in h_masses]
    ipdb.set_trace()
