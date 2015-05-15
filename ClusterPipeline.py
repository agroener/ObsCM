import numpy as np
from xlrd import open_workbook
import matplotlib.pyplot as plt
import emcee
import triangle

# imports for testing
import time
import ipdb

# For changing plotting parameters 
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

## The general conversion process will take data in one convention (either 200 or virial), and
## convert these measurements to a common cosmology and uncertainty convention.


## Functions for converting between overdensity convention (if needed) come from MConvert_Personal.py.
import MConvert_Personal as MC
## Functions for converting between cosmologies come from CosmoConvert.py.
import CosmoConvert as CC
## Functions for converting confidence intervals to symmetric uncertainties come from Uncertainties.py.
import Uncertainties as UN
## Function for modeling the normalized masses as a GMM
from GaussianMixtureModel import GMM
## Functions for loglikelihood
from LogLikelihood import lnprob

# Get the data form the excel sheet cm_data.xlsx
import DataHandler as DH
print("Importing data...")
clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology=DH.startup()

# Printing out number of total measurements, number of unique clusters, and number of studies drawn from
print("There are a total of {} measurements in the database.".format(len(clusters)))
print("There are a total of {} unique cluster objects.".format(len(set(clusters))))
print("There are a total of {} studies.".format(len(set(short_refs))))

# Take all virial measurements and normalize uncertainties
mvir_norm,mvir_p_norm,mvir_m_norm,cvir_norm,cvir_p_norm,cvir_m_norm = ([],[],[],[],[],[])
methods_norm,z_norm,cl_norm = ([],[],[])
lamb = 0.75 # the multiplicative factor which controls how much to shift the mean in the procedure for adjusting
print("Normalizing uncertainties...")
for i in range(len(clusters)):
    if mvir[i] not in [u'TBD',u'nan']:
        tmp = UN.normalize_uncertainty(mvir[i],mvir_plus[i],mvir_minus[i],cvir[i],cvir_plus[i],cvir_minus[i],lamb=lamb)
        mvir_norm.append(tmp[0])
        mvir_p_norm.append(tmp[1])
        mvir_m_norm.append(tmp[2])
        cvir_norm.append(tmp[3])
        cvir_p_norm.append(tmp[4])
        cvir_m_norm.append(tmp[5])
        methods_norm.append(methods[i])
        z_norm.append(redshift[i])
        cl_norm.append(clusters[i])

## Temporary section for outputting data for each method at this point
## Does not co-add like measurements for clusters together; for the moment
## it treats them as separate measurements.
'''
import GenDataForLinearReg as GD
print("Outputting X-ray data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='x-ray',plot=True)
print("Outputting WL data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='wl',plot=True)
print("Outputting SL data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='sl',plot=True)
print("Outputting WL+SL data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='wl+sl',plot=True)
print("Outputting CM data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='cm',plot=True)
print("Outputting LOSVD data...")
GD.writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method='losvd',plot=True)
ipdb.set_trace()
'''


## Normalize over cosmology at this point



# Find best fit parameters for GMM
weights,means,covs = GMM(mvir_norm)


# Perform parameter fitting here
'''
#lnlike([0.1,1,3],(weights,means,covs),mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,z_norm)
ndim, nwalkers = 3, 50
nsteps = 10
nbins = 100

pos = [1e-4*np.random.randn(ndim) for i in range(nwalkers)]
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                args=((weights,means,covs),
                                      mvir_norm, mvir_p_norm,
                                      cvir_norm, cvir_p_norm, z_norm),
                                threads=2)
t1 = time.time()
sampler.run_mcmc(pos, nsteps)
t2 = time.time()
print("Runtime for nwalkers={}/nsteps={}: {}".format(nwalkers,nsteps,str(t2-t1)))
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

fig = triangle.corner(samples, labels=["$\\alpha$", "$\\beta$", "$\sigma^{2}$"],
                          quantiles=[0.16, 0.5, 0.84],
                          plot_contours=True,plot_datapoints=True, bins=nbins)
plt.show()
'''



#################################################
## Some Useful Functions - Mostly For Plotting ##
#################################################
def plot_results(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,z_norm):
    # Plotting the results
    plt.figure(figsize=(8,8))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\mathrm{M_{vir}/M_{\odot}}$')
    plt.ylabel(r'$\mathrm{c_{vir}(1+z)}$')
    for i in range(len(mvir_norm)):
        if methods_norm[i] in ['WL']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='purple')
        if methods_norm[i] in ['SL']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='red')
        if methods_norm[i] in ['WL+SL']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='black')
        if methods_norm[i] in ['X-ray']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='green')
        if methods_norm[i] in ['CM']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='blue')
        if methods_norm[i] in ['LOSVD']:
            plt.errorbar(mvir_norm[i]*1e14,(1+z_norm[i])*cvir_norm[i],yerr=[cvir_p_norm[i]*(1+z_norm[i])],xerr=[mvir_p_norm[i]*1e14],color='orange')
    plt.show()

def plot_redshift_mass_distr(mvir_norm,z_norm,methods_norm):
    # Plot distribution of masses for each method
    plt.figure(figsize=(8,8))
    plt.xlabel(r'$\mathrm{M_{vir}/M_{\odot}}$',fontsize=18)
    plt.ylabel(r'$\mathrm{z}$',fontsize=18)
    plt.xscale('log')
    wl = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'WL']
    sl = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'SL']
    wl_sl = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'WL+SL']
    cm = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'CM']
    xray = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'X-ray']
    losvd = [mvir_norm[i]*1e14 for i in range(len(mvir_norm)) if methods_norm[i] == 'LOSVD']
    plt.scatter(wl,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'WL'],color='purple',marker='d',label='WL')
    plt.scatter(sl,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'SL'],color='red',marker='s',label='SL')
    plt.scatter(wl_sl,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'WL+SL'],color='black',marker='o',label='WL+SL')
    plt.scatter(cm,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'CM'],color='blue',marker='x',label='CM')
    plt.scatter(xray,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'X-ray'],color='green',marker='*',label='X-ray')
    plt.scatter(losvd,[z_norm[i] for i in range(len(z_norm)) if methods_norm[i] == 'LOSVD'],color='orange',marker='^',label='LOSVD')
    plt.legend(loc=0,scatterpoints=1,fontsize=10)
    plt.xlim(1e13,3e16)
    plt.ylim(0,1.5)
    plt.show()

def compare_methods(method1, method2, scaleaxes = True, scaleto = None):

    # Setting up method 1 data
    m1_cl = [cl_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_concs = [cvir_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_concs_err = [cvir_p_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_mass = [mvir_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_mass_err = [mvir_p_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_z = [z_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method1]
    m1_cl_n, m1_concs_n, m1_concs_err_n, m1_mass_n, m1_mass_err_n, m1_z_n = ([],[],[],[],[],[])
    for i in set(m1_cl):
        # handling co-adding measurements together from repeat clusters
        if m1_cl.count(i) > 1:
            tmp_indices = [j for j in range(len(m1_cl)) if m1_cl[j] == i]
            m_rep = [m1_mass[j] for j in tmp_indices]
            m_err_rep = [m1_mass_err[j] for j in tmp_indices]
            m_weights = [(1.0/m_err_rep[j]**2) for j in range(len(m_rep))]
            m_new = sum([m_rep[j]*m_weights[j] for j in range(len(m_rep))])/sum(m_weights)
            m_err_new = 1.0/np.sqrt(sum(m_weights))
            c_rep = [m1_concs[j] for j in tmp_indices]
            c_err_rep = [m1_concs_err[j] for j in tmp_indices]
            c_weights = [(1.0/c_err_rep[j]**2) for j in range(len(c_rep))]
            c_new = sum([c_rep[j]*c_weights[j] for j in range(len(c_rep))])/sum(c_weights)
            c_err_new = 1.0/np.sqrt(sum(c_weights))
            m1_cl_n.append(i)
            m1_concs_n.append(c_new)
            m1_concs_err_n.append(c_err_new)
            m1_mass_n.append(m_new)
            m1_mass_err_n.append(m_err_new)
            m1_z_n.append(m1_z[m1_cl.index(i)])
        else:
            tmp_index = m1_cl.index(i)
            m1_cl_n.append(i)
            m1_concs_n.append(m1_concs[tmp_index])
            m1_concs_err_n.append(m1_concs_err[tmp_index])
            m1_mass_n.append(m1_mass[tmp_index])
            m1_mass_err_n.append(m1_mass_err[tmp_index])
            m1_z_n.append(m1_z[tmp_index])

    # Setting up method 2 data
    m2_cl = [cl_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_concs = [cvir_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_concs_err = [cvir_p_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_mass = [mvir_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_mass_err = [mvir_p_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_z = [z_norm[i] for i in range(len(cl_norm)) if methods_norm[i] == method2]
    m2_cl_n, m2_concs_n, m2_concs_err_n, m2_mass_n, m2_mass_err_n, m2_z_n = ([],[],[],[],[],[])
    for i in set(m2_cl):
        # handling co-adding measurements together from repeat clusters
        if m2_cl.count(i) > 1:
            tmp_indices = [j for j in range(len(m2_cl)) if m2_cl[j] == i]
            m_rep = [m2_mass[j] for j in tmp_indices]
            m_err_rep = [m2_mass_err[j] for j in tmp_indices]
            m_weights = [(1.0/m_err_rep[j]**2) for j in range(len(m_rep))]
            m_new = sum([m_rep[j]*m_weights[j] for j in range(len(m_rep))])/sum(m_weights)
            m_err_new = 1.0/np.sqrt(sum(m_weights))
            c_rep = [m2_concs[j] for j in tmp_indices]
            c_err_rep = [m2_concs_err[j] for j in tmp_indices]
            c_weights = [(1.0/c_err_rep[j]**2) for j in range(len(c_rep))]
            c_new = sum([c_rep[j]*c_weights[j] for j in range(len(c_rep))])/sum(c_weights)
            c_err_new = 1.0/np.sqrt(sum(c_weights))
            m2_cl_n.append(i)
            m2_concs_n.append(c_new)
            m2_concs_err_n.append(c_err_new)
            m2_mass_n.append(m_new)
            m2_mass_err_n.append(m_err_new)
            m2_z_n.append(m2_z[m2_cl.index(i)])
        else:
            tmp_index = m2_cl.index(i)
            m2_cl_n.append(i)
            m2_concs_n.append(m2_concs[tmp_index])
            m2_concs_err_n.append(m2_concs_err[tmp_index])
            m2_mass_n.append(m2_mass[tmp_index])
            m2_mass_err_n.append(m2_mass_err[tmp_index])
            m2_z_n.append(m2_z[tmp_index])

    # Finding repeat clusters with measurements which appear in method 1 and method 2
    both_cl = [i for i in set(m1_cl).intersection(m2_cl)]
    both_z1 = [m1_z_n[i] for i in range(len(m1_cl_n)) if m1_cl_n[i] in both_cl]
    both_z2 = [m2_z_n[i] for i in range(len(m2_cl_n)) if m2_cl_n[i] in both_cl]
    both_z = [i for i in set(both_z1 + both_z2)]
    
    # Plotting the results
    plt.figure(figsize=(7,7))
    plt.title("({}/{})".format(method1,method2),fontsize=14)
    plt.xlabel(r'$\mathrm{\frac{M_{vir}}{M_{vir}} - 1}$',fontsize=20)
    plt.ylabel(r'$\mathrm{\frac{c_{vir}}{c_{vir}} - 1}$',fontsize=20,rotation='horizontal')
    if scaleaxes is True:
        if scaleto is not None:
            plt.xlim(-1,scaleto)
            plt.ylim(-1,scaleto)
        else:
            plt.xlim(-1,5)
            plt.ylim(-1,5)
        filename = "{}_{}_Comparison_Scaled.png".format(method1, method2)
    else:
        filename = "{}_{}_Comparison_NotScaled.png".format(method1, method2)
    plt.axvline(x=0,color='black',linestyle='--')
    plt.axhline(y=0,color='black',linestyle='--')
    maxz = max(both_z)
    cm = plt.cm.get_cmap('RdYlBu')
    for it, vals in enumerate(both_cl):
        m1_index = m1_cl_n.index(vals)
        m2_index = m2_cl_n.index(vals)
        tmp_x = m1_mass_n[m1_index]/m2_mass_n[m2_index]
        tmp_y = m1_concs_n[m1_index]/m2_concs_n[m2_index]
        tmp_sigx = abs(m1_mass_n[m1_index]/m2_mass_n[m2_index]-1) * (m1_mass_err_n[m1_index]/m1_mass_n[m1_index] + m2_mass_err_n[m2_index]/m2_mass_n[m2_index])
        tmp_sigy = abs(m1_concs_n[m1_index]/m2_concs_n[m2_index]-1) * (m1_concs_err_n[m1_index]/m1_concs_n[m1_index] + m2_concs_err_n[m2_index]/m2_concs_n[m2_index])
        plt.errorbar(tmp_x-1,tmp_y-1,yerr=tmp_sigy,xerr=tmp_sigx, zorder=1, color='black')
        sc = plt.scatter(tmp_x-1, tmp_y-1, c=m1_z_n[m1_index], vmin=0, vmax=maxz, s=43, cmap=cm, zorder=2)
    plt.colorbar(sc)
    # Saves figure onto Desktop; path should work for laptop and work computer only
    plt.savefig("/Users/groenera/Desktop/{}".format(filename))
    plt.show()

    
            
if __name__ == "__main__":
    
    print()
    #plot_redshift_mass_distr(mvir_norm,z_norm,methods_norm)

    ## Plotting concentration/mass comparisons of clusters which are measured in both
    #compare_methods('X-ray', 'WL',scaleaxes=True,scaleto=5)
    #compare_methods('WL', 'WL+SL',scaleaxes=True,scaleto=3)
    #compare_methods('CM', 'LOSVD',scaleaxes=True,scaleto=5)
