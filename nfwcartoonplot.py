import numpy as np
import matplotlib.pyplot as plt

import ipdb

def nfw(r,c,rvir=1.0):
    rs = rvir/c
    delta_c = (200/3)*(c**3/(np.log(1+c)-(c/(1+c))))
    return delta_c/((r/rs)*(1+(r/rs)))



if __name__ == "__main__":
    f, ax = plt.subplots()
    rlist = np.linspace(0.01,20,1000)
    #ipdb.set_trace()
    ax.plot(np.log10(rlist),np.log10(nfw(rlist,2.0)),label=r'$\mathrm{c=2}$',color='blue')
    ax.axvline(x=np.log10(1.0/2.0),ymax=0.54,linewidth=3,linestyle='--',color='blue')
    ax.plot(np.log10(rlist),np.log10(nfw(rlist,4.0)),label=r'$\mathrm{c=4}$',color='green')
    #ax.axvline(x=np.log10(1.0/4.0),ymax=0.63,linewidth=3,linestyle='--',color='green')
    ax.plot(np.log10(rlist),np.log10(nfw(rlist,8.0)),label=r'$\mathrm{c=8}$',color='red')
    #ax.axvline(x=np.log10(1.0/8.0),ymax=0.73,linewidth=3,linestyle='--',color='red')
    ax.plot(np.log10(rlist),np.log10(nfw(rlist,16.0)),label=r'$\mathrm{c=16}$',color='cyan')
    #ax.axvline(x=np.log10(1.0/16.0),ymax=0.845,linewidth=3,linestyle='--',color='cyan')
    ax.set_ylabel(r"$\mathrm{\log\,\frac{\rho}{\rho_{cr}}}$",fontsize=20)
    ax.set_xlabel(r"$\mathrm{\log(r)}$",fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    labels.reverse()
    handles.reverse()
    ax.legend(handles,labels,loc=0)
    plt.tight_layout()
    plt.show()
