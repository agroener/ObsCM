from __future__ import division # float division by default
import numpy as np
from scipy.integrate import quad
import ipdb
import MConvert_Personal as mc
import matplotlib.pyplot as plt

def MConvert_SE14(Mold,DeltaOld,DeltaNew,ConcentOld,redshift,method,Omega_m_0_new,Omega_L_0_new,Omega_m_0_old=0.3,Omega_L_0_old=0.7):
    Mnew = mc.Mconvert(Mold,DeltaOld,DeltaNew,ConcentOld)
    Cnew = mc.Cconvert(Mold,DeltaOld,DeltaNew,ConcentOld)
    if method is None:
        return Mnew,Cnew
    elif method in ['WL','SL','WL+SL']:
        z_source = 1
        # original cosmology
        Dds_c1 = DAng(1/(1+z_source),1/(1+redshift),Omega_m_0_old,Omega_L_0_old)
        Ds_c1 = DAng(1/(1+z_source),1.0,Omega_m_0_old,Omega_L_0_old)
        H_c1 = Hubble(redshift,Omega_m_0_old,Omega_L_0_old)
        # new cosmology
        Dds_c2 = DAng(1/(1+z_source),1/(1+redshift),Omega_m_0_new,Omega_L_0_new)
        Ds_c2 = DAng(1/(1+z_source),1.0,Omega_m_0_new,Omega_L_0_new)
        H_c2 = Hubble(redshift,Omega_m_0_new,Omega_L_0_new)
        # conversion ratio
        R = (((Dds_c2/Ds_c2)**(-3.0/2))/((Dds_c1/Ds_c1)**(-3.0/2)))*(H_c2/H_c1)
        return R

def conversion_ratio(redshift,method,Omega_m_0_new,Omega_L_0_new,Omega_m_0_old=0.3,Omega_L_0_old=0.7):
    z_source = 1
    # original cosmology
    Dds_c1 = DAng(1/(1+z_source),1/(1+redshift),Omega_m_0_old,Omega_L_0_old)
    Ds_c1 = DAng(1/(1+z_source),1.0,Omega_m_0_old,Omega_L_0_old)
    H_c1 = Hubble(redshift,Omega_m_0_old,Omega_L_0_old)
    # new cosmology
    Dds_c2 = DAng(1/(1+z_source),1/(1+redshift),Omega_m_0_new,Omega_L_0_new)
    Ds_c2 = DAng(1/(1+z_source),1.0,Omega_m_0_new,Omega_L_0_new)
    H_c2 = Hubble(redshift,Omega_m_0_new,Omega_L_0_new)
    # conversion ratio
    R = (((Dds_c2/Ds_c2)**(-3.0/2))/((Dds_c1/Ds_c1)**(-3.0/2)))*(H_c2/H_c1)
    return R
    
def DAng(a_source,a_cluster=1,Omega_M=0.3,Omega_L=0.7): # For flat cosmologies only 
    return Chi(a_source,a_cluster,Omega_M,Omega_L)*a_source # in units of c/H0

def Chi(a_source,a_cluster=1,Omega_M=0.3,Omega_L=0.7):
    I = quad(lambda A: 1.0/((A**2)*np.sqrt(Omega_M/(A**3) + Omega_L)),a_source,a_cluster)
    return I[0] # in units of c/H0
    
def Hubble(z,Omega_M=0.3,Omega_L=0.7,Omega_k=0.0,Omega_g=0.0,h=0.7):
    return 100*h*np.sqrt(Omega_M*((1+z)**3) + Omega_L + Omega_k*((1+z)**2) + Omega_g*((1+z)**4)) # in km/s/Mpc


if __name__ == "__main__":
    #MConvert_SE14(10.5,200,500,7.5,0.183,'WL',0.32,0.68)

    # Making the M(c1)/M(c2) plot
    omega_m_list = np.linspace(0.0,1.0,1000)
    rat = [conversion_ratio(0.183,'WL',i,1-i) for i in omega_m_list]
    rat_u = [abs(rat[i]-1.01) for i in range(len(rat))]
    lim_u = rat_u.index(min(rat_u))
    rat_l = [abs(rat[i]-0.99) for i in range(len(rat))]
    lim_l = rat_l.index(min(rat_l))
    plt.plot(omega_m_list,rat)
    plt.title(r'$\Omega_{m,0}=0.3$, $\Omega_{\Lambda,0}=0.7$, $h=0.7$')
    plt.xlabel(r'$\Omega_{m,0}$'+" " + r'$(\Omega_{\Lambda,0} = 1 - \Omega_{m,0})$',fontsize=16)
    plt.ylabel(r'$ M_{\Delta}(c1)/M_{\Delta}(c2) = \left( \frac{D_{ds}(c2)/D_{s}(c2)}{D_{ds}(c1)/D_{s}(c1)} \right)^{-3/2} H(c2)/H(c1)$',fontsize=17)
    plt.axvline(x=0.3,color='red',linewidth=3,linestyle="--")
    plt.axvspan(omega_m_list[lim_l],omega_m_list[lim_u],alpha=0.5,color='red')
    plt.axhspan(rat[lim_l],rat[lim_u],alpha=0.5,color='blue')
    plt.text(0.1,1.02,"1.01")
    plt.text(0.1,0.97,"0.99")
    plt.text(0.12,0.85,r"$\Omega_{m,0} = $"+"{}".format(round(omega_m_list[lim_l],3),rotation='horizontal'))
    plt.text(0.33,0.85,r"$\Omega_{m,0} = $"+"{}".format(round(omega_m_list[lim_u],3),rotation='horizontal'))
    plt.show()
