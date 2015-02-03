# For the data
from xlrd import open_workbook

#For MCMC fitting
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as op
import emcee
import triangle
from matplotlib.colors import LogNorm

# Debugging
import ipdb

# Opening Excel File
#wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx') # for laptop and work machine
wb = open_workbook('/home/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx') # for home desktop

# First sheet
sheet_names = wb.sheet_names()
sheet1 = wb.sheet_by_name(sheet_names[0])

# Get headers and data
headers = []
clusters = sheet1.col_values(0)
headers.append(clusters.pop(0))
redshift = sheet1.col_values(1)
headers.append(redshift.pop(0))
methods = sheet1.col_values(2)
headers.append(methods.pop(0))
c200 = sheet1.col_values(3)
headers.append(c200.pop(0))
c200_plus = sheet1.col_values(4)
headers.append(c200_plus.pop(0))
c200_minus = sheet1.col_values(5)
headers.append(c200_minus.pop(0))
m200 = sheet1.col_values(6)
headers.append(m200.pop(0))
m200_plus = sheet1.col_values(7)
headers.append(m200_plus.pop(0))
m200_minus = sheet1.col_values(8)
headers.append(m200_minus.pop(0))
cvir = sheet1.col_values(9)
headers.append(cvir.pop(0))
cvir_plus = sheet1.col_values(10)
headers.append(cvir_plus.pop(0))
cvir_minus = sheet1.col_values(11)
headers.append(cvir_minus.pop(0))
mvir = sheet1.col_values(12)
headers.append(mvir.pop(0))
mvir_plus = sheet1.col_values(13)
headers.append(mvir_plus.pop(0))
mvir_minus = sheet1.col_values(14)
headers.append(mvir_minus.pop(0))
short_refs = sheet1.col_values(15)
headers.append(short_refs.pop(0))
orig_convention = sheet1.col_values(16)
headers.append(orig_convention.pop(0))
cosmology = sheet1.col_values(17)
headers.append(cosmology.pop(0))

mvir_tbd = [i for i in range(len(mvir)) if mvir[i] == 'TBD' and methods[i] == 'CM']
cvir_tbd = [i for i in range(len(cvir)) if cvir[i] == 'TBD' and methods[i] == 'CM']
m200_tbd = [i for i in range(len(m200)) if m200[i] == 'TBD' and methods[i] == 'CM']
c200_tbd = [i for i in range(len(c200)) if c200[i] == 'TBD' and methods[i] == 'CM']

assert len(mvir) == len(cvir)
assert mvir_tbd == cvir_tbd 
assert len(m200) == len(c200)
assert m200_tbd == c200_tbd
assert len(mvir) == len(m200)
assert mvir_tbd == m200_tbd

mvir_cm = [mvir[i]*1e14 for i in range(len(mvir)) if i not in mvir_tbd and methods[i] == 'CM']
cvir_cm = [cvir[i] for i in range(len(cvir)) if i not in cvir_tbd and methods[i] == 'CM']
m200_cm = [m200[i]*1e14 for i in range(len(m200)) if i not in m200_tbd and methods[i] == 'CM']
c200_cm = [c200[i] for i in range(len(c200)) if i not in c200_tbd and methods[i] == 'CM']

z_cm = [redshift[i] for i in range(len(redshift)) if i not in mvir_tbd and methods[i] == 'CM']

#___________________________________________
# Functions #
############################################
def lnlike(theta, x, y, z):
    alpha, c0 = theta
    model = ((x/mstar)**alpha) * (c0/(1+z))
    return -0.5*(np.sum((y-model)**2))

def lnprior(theta):
    alpha, c0 = theta
    if 0.0 < alpha < 0.5 and 0.0 < c0 < 20.0:
        return 0.0
    return -np.inf

def lnprob(theta, x, y, z):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, z)

nll = lambda *args: -lnlike(*args)
#___________________________________________

# constants
mstar = 1.857e13 #1.3e13Msun/h; where h=07
ndim, nwalkers = 2, 1000
nbins = 100

# initializing walker positions
pos = [1e-4*np.random.randn(ndim) for i in range(nwalkers)]

# setting up and running the algorithm
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, a=2.0, 
                                args=(np.array(mvir_cm), np.array(cvir_cm), np.array(z_cm)),
                                threads=2)
sampler.run_mcmc(pos, 750)
samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

if __name__ == '__main__':

    # plotting the results
    fig = triangle.corner(samples, labels=["$\\alpha$", "$C_{0}$"],
                          quantiles=[0.16, 0.5, 0.84],
                          plot_contours=True,plot_datapoints=True, bins=nbins)
    plt.show()

    # plotting the fit over the data
    xl = np.array([min(mvir_cm)*1.10, max(mvir_cm)*1.10])
    for alpha, c0 in samples[np.random.randint(len(samples), size=300)]:
        plt.plot(xl, ((xl/mstar)**alpha) * (c0/(1+0.01)), color="k", alpha=0.1)
    plt.scatter(mvir_cm, cvir_cm)
    plt.show()
