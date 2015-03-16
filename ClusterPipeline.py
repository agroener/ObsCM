import numpy as np
from xlrd import open_workbook
import matplotlib.pyplot as plt
import emcee
import triangle

# imports for testing
import time
import ipdb

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

### Normalize over cosmology at this point
        
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

# Find best fit parameters for GMM
weights,means,covs = GMM(mvir_norm)

# Perform parameter fitting here

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

