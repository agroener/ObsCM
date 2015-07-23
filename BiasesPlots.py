from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq

import ipdb

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


if __name__ == "__main__":

    ## Plotting
    # C/c plot
    plt.figure(1)
    conclist = [1,3,5,7,9]
    shapelist = np.linspace(0.25,1.0,100)
    outlist = [[(conc_finder_pro(conc,[0.0],shape)[0][0])/conc for shape in shapelist] for conc in conclist]
    plt.ylabel(r"$\mathrm{\frac{C_{2D}}{c}}$",fontsize=18,rotation="horizontal")
    plt.xlabel(r"$q$",fontsize=18)
    for i in range(len(outlist)):
        plt.plot(shapelist,outlist[i],label="c={}".format(conclist[i]))
    plt.xlim(0.25,1.0)
    plt.legend(loc=0)
    # Q vs q and theta plot
    plt.figure(2)
    qlist = [0.25,0.5,0.75,1.0]
    thetalist = np.linspace(0.0,np.pi,100)
    outlist2 = [[Q(shape,shape,0.0,theta) for theta in thetalist] for shape in qlist]
    plt.ylabel(r"$\mathrm{Q}$",fontsize=18,rotation="horizontal")
    plt.xlabel(r"$\mathrm{\theta}$",fontsize=18)
    for i in range(len(outlist2)):
        plt.plot(thetalist,outlist2[i],label="q={}".format(qlist[i]))
    plt.legend(loc=0)
    plt.ylim(-0.1,1.1)
    plt.xlim(0,np.pi)
    plt.axvline(x=np.pi/2,linestyle='--',color='black',linewidth=2)
    plt.text(np.pi/2,0.6,"Perpendicular to Major Axis",rotation=90,fontsize=17)
    plt.show()
