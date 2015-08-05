# A script to generate something like Figure 12 in Prada et al. (2012)
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

import MConvert_Personal as MC

import ipdb

Omega_m_0 = 0
Omega_L_0 = 0

def set_cosmology_global(Omega_matter_0=0.3,Omega_Lambda_0=0.7):
    global Omega_m_0
    global Omega_L_0
    Omega_m_0 = Omega_matter_0
    Omega_L_0 = Omega_Lambda_0

def D(a,Omega_m_0,Omega_L_0):
    #print(Omega_m_0,Omega_L_0)
    x = (Omega_L_0/Omega_m_0)**(1/3) * a
    func = lambda t: (t**(3/2))/((1+t**3)**(3/2))
    Integral = quad(func, 0, x)
    return (5/2) * (Omega_m_0/Omega_L_0)**(1/3) * ( (np.sqrt(1+x**3)) / (x**(3/2))) * (Integral[0])

def c(M,z):
    a = 1/(1+z)
    x = (0.7/0.3)**(1/3) * a
    sigma_p = sigma_prime(M,a)
    return B0(x)*C(sigma_p)

def B0(x):
    return c_min(x)/c_min(1.393)

# Compared against Fig. 10 in paper. - Looks good.
def c_min(x):
    return 3.681 + (5.033 - 3.681) * ( (1/np.pi) * np.arctan(6.948*(x - 0.424)) + (1/2) )

def B1(x):
    return sigma_inv_min(x)/sigma_inv_min(1.393)

def C(sigma_prime):
    return 2.881 * ( (sigma_prime/1.257)**1.022 + 1 ) * (np.exp(0.060/(sigma_prime**2)))

def sigma_prime(M,a):
    x = (0.7/0.3)**(1/3) * a
    return B1(x)*sigma(M,a) 

def sigma(M,a):
    y = 1/(M*0.7/10**12)
    return D(a,Omega_m_0,Omega_L_0) * ( (16.9*(y**0.41)) / (1 + 1.102*(y**0.20) + 6.22*(y**0.333)) ) 

# Compared against Fig. 10 in paper. - Looks good.
def sigma_inv_min(x):
    return 1.047 + (1.646 - 1.047) * ( (1/np.pi) * np.arctan(7.386*(x-0.526)) + (1/2) )

def PradaRelation(Omega_m_0=0.3,Omega_L_0=0.7,z=0.0):
    # Set cosmology
    set_cosmology_global(Omega_m_0,Omega_L_0)
    # Create mass list; calculate concentrations
    M = np.linspace(10**11,5*10**15,200)
    cout = c(M,z)
    # Convert values
    #ipdb.set_trace()
    cvir = MC.Cconvert(M,200,MC.DeltaFinder(Omega_m_0,Omega_L_0,z),cout)
    mvir = MC.Mconvert(M,200,MC.DeltaFinder(Omega_m_0,Omega_L_0,z),cout)
    return cvir,mvir
    
    
if __name__ == "__main__":
    # Producing the upturn plot
    '''
    set_cosmology_global(0.30,0.70)
    print(Omega_m_0,Omega_L_0)
    M1 = np.linspace(10**11,10**15,200)
    cout1 = c(M1,z=0.0)
    M2 = np.linspace(10**11,10**14.6,200)
    cout2 = c(M2,z=0.5)
    M3 = np.linspace(10**11,10**14.4,200)
    cout3 = c(M3,z=1.0)
    M4 = np.linspace(10**11,10**14,200)
    cout4 = c(M4,z=2.0)
    f, axarr = plt.subplots(2, sharex=True)
    axarr[0].plot(np.log10(M1),np.log10(cout1),linewidth=2,color='black',label='z=0.0',linestyle=':')
    axarr[0].plot(np.log10(M2),np.log10(cout2),linewidth=2,color='yellow',label='z=0.5')
    axarr[0].plot(np.log10(M3),np.log10(cout3),linewidth=2,color='green',label='z=1.0',linestyle='-.')
    axarr[0].plot(np.log10(M4),np.log10(cout4),linewidth=2,color='red',label='z=2.0',linestyle='--')
    axarr[0].set_ylabel(r'$\mathrm{\log \, c}$',fontsize=20)
    axarr[1].plot(np.log10(M1),np.log10(cout1*(1+0.0)),linewidth=2,color='black',label='z=0.0',linestyle=':')
    axarr[1].plot(np.log10(M2),np.log10(cout2*(1+0.5)),linewidth=2,color='yellow',label='z=0.5')
    axarr[1].plot(np.log10(M3),np.log10(cout3*(1+1.0)),linewidth=2,color='green',label='z=1.0',linestyle='-.')
    axarr[1].plot(np.log10(M4),np.log10(cout4*(1+2.0)),linewidth=2,color='red',label='z=2.0',linestyle='--')
    axarr[1].set_xlabel(r'$\mathrm{\log \, M \, (h^{-1} M_{\odot})}$',fontsize=20)
    axarr[1].set_ylabel(r'$\mathrm{\log \, c \, (1+z)}$',fontsize=20)
    axarr[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
          ncol=2, fancybox=True, shadow=True)
    plt.tight_layout()
    plt.show()
    '''
    # Relation for comparing with WL and WL+SL relations
    PradaRelation(0.3,0.7,0.5)
