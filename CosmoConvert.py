from __future__ import division
import matplotlib.pyplot as plt
from scipy.integrate import *
import numpy as np
import ipdb

# For changing plotting parameters 
from matplotlib import rcParams
#rcParams.update({'figure.autolayout': True})

#######################
### Constants (mks) ###
#######################

G = 6.673e-11 # m^3 kg^-1 s^-2
c = 3e8 

#################
### Functions ###
#################

# all functions assume a flat cosmology (Omega_M + Omega_L = 1)

def Hubble(z,Omega_M=0.3,Omega_L=0.7,h=0.7):
    H0 = 100*h*3.24077929e-20 #mks
    return H0*np.sqrt(Omega_M*((1+z)**3) + Omega_L) 

def ComovingDistance(z,Omega_M=0.3,Omega_L=0.7):
    H0 = Hubble(0,Omega_M,Omega_L)
    dH = 3e8/H0
    I = quad(lambda Z: 1/(np.sqrt(Omega_M*((1+Z)**3) + Omega_L)),0,z)
    return I[0]*dH

def AngularDiameterDistance(z_low,z_high,Omega_M=0.3,Omega_L=0.7):
    Chi_low = ComovingDistance(z_low,Omega_M,Omega_L)
    Chi_high = ComovingDistance(z_high,Omega_M,Omega_L)
    return (1/(1+z_high)) * (Chi_high - Chi_low)

def SigmaCrit(z_lens,z_source,Omega_M=0.3,Omega_L=0.7):
    if z_source == 'inf': # source at infinity
        return 3e8 * 3e8 / (4 * np.pi * G * AngularDiameterDistance(0,z_lens,Omega_M=Omega_M,Omega_L=Omega_L) )
    else:
        return c**2 * AngularDiameterDistance(0,z_source,Omega_M,Omega_L) / (4 * np.pi * G * AngularDiameterDistance(0,z_lens,Omega_M,Omega_L) * AngularDiameterDistance(z_lens,z_source,Omega_M,Omega_L) )

def T(z_lens,z_source,OM1,OL1,OM2,OL2):
    return SigmaCrit(z_lens,z_source,OM2,OL2)/SigmaCrit(z_lens,z_source,OM1,OL1)

def R(OM1,OL1,Delta1,OM2,OL2,Delta2,z_lens,h=0.7):
    num = AngularDiameterDistance(0,z_lens,OM1,OL1)**3 * Hubble(z,OM1,OL1,h)**2 * Delta1
    denom = AngularDiameterDistance(0,z_lens,OM2,OL2)**3 * Hubble(z,OM2,OL2,h)**2 * Delta2
    return num/denom

def f(x):
    return np.log(1+x) - x/(1+x)

# lensing only 
def normalize_cosmology(m,m_p,m_m,c,c_p,c_m,OM1,OL1,OM2,OL2,z):
    RHS = (f(c200)*T(z,'inf',OM1,OL1,OM2,OL2)) / (c200**3*R(OM1,OL1,Delta1,OM2,OL2,Delta2,z))
    
    

if __name__ == "__main__":
    # Exploring the effect of choosing source redshift on conversion factor calculations
    '''
    z=0.1
    h=0.7
    zsource = np.linspace(0.2,100,1000)
    list1 = [AngularDiameterDistance(z,zsource[i],0.3,0.7) for i in range(len(zsource))]
    list2 = [AngularDiameterDistance(z,zsource[i],0.27,0.73) for i in range(len(zsource))]
    list3 = [AngularDiameterDistance(0,zsource[i],0.3,0.7) for i in range(len(zsource))]
    list4 = [AngularDiameterDistance(0,zsource[i],0.27,0.73) for i in range(len(zsource))]
    list5 = [T(z,i,0.3,0.7,0.27,0.73) for i in zsource]
    list6 = [T(z,i,0.3,0.7,0.25,0.75) for i in zsource]
    list7 = [T(z,i,0.3,0.7,0.33,0.67) for i in zsource]
    difflist1 = [(list1[i]/list2[i])*100 for i in range(len(zsource))]
    difflist2 = [(list3[i]/list4[i])*100 for i in range(len(zsource))]
    plt.figure(figsize=(8,8))
    plt.plot(zsource,difflist1,label=r'$\mathrm{D_{ds}(z_{L}=0.1,z_{S}=z)}$')
    plt.plot(zsource,difflist2,label=r'$\mathrm{D_{s}(z_{S}=z)}$')
    plt.xlabel(r'$\mathrm{z}$',fontsize=18)
    plt.ylabel(r'$\mathrm{\left( \frac{D_{A}\Omega_{1}}{D_{A}\Omega_{2}} \right)\times 100}$',fontsize=18)
    plt.legend(loc=0)
    plt.show()
    plt.figure(figsize=(8,8))
    plt.plot(zsource,list5,label=r'$\Omega_{2}=\{ \Omega_{m}=0.27,\Omega_{\Lambda}=0.73 \}$')
    plt.plot(zsource,list6,label=r'$\Omega_{2}=\{ \Omega_{m}=0.25,\Omega_{\Lambda}=0.75 \}$')
    plt.plot(zsource,list7,label=r'$\Omega_{2}=\{ \Omega_{m}=0.33,\Omega_{\Lambda}=0.67 \}$')
    plt.xlabel(r'$\mathrm{z}$',fontsize=18)
    plt.ylabel(r'$T = \frac{\Sigma_{cr} \Omega_{2}}{\Sigma_{cr} \Omega_{1}}$',fontsize=18)
    plt.legend(loc=0)
    plt.show()
    '''

    # Calculating RHS of equation for example in pdf
    '''
    z = 0.1
    z_source = 20
    c200 = 2.5
    M200 = 5.1e14
    OM2 = 0.3
    OL2 = 0.7
    OM1 = 1.00
    OL1 = 0.00
    Delta1 = 200
    Delta2 = 200
    
    RHS = (f(c200)*T(z,'inf',OM1,OL1,OM2,OL2)) / (c200**3*R(OM1,OL1,Delta1,OM2,OL2,Delta2,z))

    clist = np.linspace(0.1,20,10000)
    LHS = [abs((f(i)/i**3)-RHS) for i in clist]
    c200_new = clist[LHS.index(min(LHS))]
    M200_new = (f(c200_new)*M200) / (f(c200)*T(z,'inf',OM1,OL1,OM2,OL2))
    '''

    # Plotting the cosmology correction for a range of cosmologies
    z = 0.1
    z_source = 20
    # testing this out with different concentrations/masses
    c200 = [1.0,5.0,9.0]
    M200 = 1.0e14
    # fiducial cosmology
    OM1 = 0.3
    OL1 = 0.7
    # new cosmologies
    OM2 = np.linspace(0.0,1.0,40)
    OL2 = [1.0-OM2[i] for i in range(len(OM2))]
    # same definition of overdensity
    Delta1 = 200
    Delta2 = 200
    # doing concentrations first
    #'''
    fig, axarr = plt.subplots(2, sharex=True)
    fig.subplots_adjust(hspace=0)
    fig.set_dpi(120)
    m200_static = 1.0e14
    for i,val in enumerate(c200):
        c200_new_list = []
        print("Iteration: {}".format(i+1))
        for j in range(len(OM2)):
            RHS = (f(val)*T(z,z_source,OM1,OL1,OM2[j],OL2[j])) / (val**3*R(OM1,OL1,Delta1,OM2[j],OL2[j],Delta2,z))
            clist = np.linspace(0.1,20,1e5)
            LHS = [abs((f(i)/i**3)-RHS) for i in clist]
            c200_new = clist[LHS.index(min(LHS))]
            c200_new_list.append(c200_new)
        c200_new_rat = [c200_new_list[i]/val for i in range(len(c200_new_list))]
        axarr[0].plot(OM2,c200_new_rat,label=r"$\mathrm{c_{200}}$"+" = {}".format(int(val)))
    axarr[0].set_ylabel(r"$\mathrm{\frac{c_{200}(\Omega_{m})}{c_{200}(\Omega_{m}=0.3)}}$",fontsize=30)
    axarr[0].axvline(x=0.3,color='black',linestyle='--',linewidth=2)
    axarr[0].axhline(y=1.0,color='black',linestyle='--',linewidth=2)
    #'''
    # doing mass second (but it's still plotted for a few different values of concentrations)
    #'''
    for i,val in enumerate(c200):
        M200_new_list = []
        print("Iteration: {}".format(i+1))
        for j in range(len(OM2)):
            RHS = (f(val)*T(z,z_source,OM1,OL1,OM2[j],OL2[j])) / (val**3*R(OM1,OL1,Delta1,OM2[j],OL2[j],Delta2,z))
            clist = np.linspace(0.1,20,1e5)
            LHS = [abs((f(k)/k**3)-RHS) for k in clist]
            c200_new = clist[LHS.index(min(LHS))]
            M200_new = (f(c200_new)*M200) / (f(val)*T(z,z_source,OM1,OL1,OM2[j],OL2[j]))
            M200_new_list.append(M200_new)
        M200_new_rat = [M200_new_list[i]/M200 for i in range(len(M200_new_list))]
        axarr[1].plot(OM2,M200_new_rat,label=r"$\mathrm{c_{200}}$"+" = {}".format(val))
    axarr[1].set_ylabel(r"$\mathrm{\frac{M_{200}(\Omega_{m})}{M_{200}(\Omega_{m}=0.3)}}$",fontsize=30)
    axarr[1].set_xlabel(r"$\mathrm{\Omega_{m} = 1 - \Omega_{\Lambda}}$",fontsize=27)
    axarr[1].axvline(x=0.3,color='black',linestyle='--',linewidth=2)
    axarr[1].axhline(y=1.0,color='black',linestyle='--',linewidth=2)
    axarr[1].legend(loc=0)
    plt.tight_layout()
    plt.show()
    #'''
