import numpy as np
import matplotlib.pyplot as plt

import ipdb

def nfw(r,c,rvir=1.0):
    rs = rvir/c
    delta_c = (200/3)*(c**3/(np.log(1+c)-(c/(1+c))))
    return delta_c/((r/rs)*(1+(r/rs)))



if __name__ == "__main__":
    rlist = np.linspace(0.01,20,1000)
    plt.plot(np.log10(rlist),np.log10(nfw(rlist,2.0)),label=r'$\mathrm{c=2}$')
    plt.plot(np.log10(rlist),np.log10(nfw(rlist,4.0)),label=r'$\mathrm{c=4}$')
    plt.plot(np.log10(rlist),np.log10(nfw(rlist,8.0)),label=r'$\mathrm{c=8}$')
    plt.plot(np.log10(rlist),np.log10(nfw(rlist,16.0)),label=r'$\mathrm{c=16}$')
    plt.ylabel(r"$\mathrm{\rho/\rho_{cr}}$",rotation='horizontal',fontsize=20)
    plt.xlabel(r"$\mathrm{r}$",fontsize=20)
    plt.legend(loc=0)
    plt.show()
    ipdb.set_trace()
