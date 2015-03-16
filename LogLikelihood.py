import numpy as np
import time

# temporary imports
import ipdb

## 50 walkers/ 10 steps = 125.3 seconds
## 
def lnlike(theta, psi, Mvir, Mvir_err, cvir, cvir_err, redshift):
    t1 = time.time()
    alpha, beta, var = theta
    weights, means, covs = psi
    weights = weights.reshape(len(weights),1)
    assert np.shape(weights) == np.shape(means)
    assert np.shape(weights) == np.shape(covs)
    x = np.array(Mvir)
    assert len(cvir) == len(redshift)
    y = np.array([cvir[i]*(1+redshift[i]) for i in range(len(cvir))])

    n = len(Mvir)
    k = len(weights)

    sum = 0
    for i in range(n):
        for j in range(k):
            V = np.matrix([[beta**2 * covs[j][0] + var + ((1+redshift[i])*cvir_err[i])**2, beta * covs[j][0]], [beta * covs[j][0], covs[j][0]+Mvir_err[i]]])
            z = np.array([[(1+redshift[i])*cvir[i]],[Mvir[i]]])
            zeta = np.array([[alpha+beta*means[j][0]],[means[j][0]]])
            sum += (weights[j][0]/(2*np.pi*np.sqrt(np.linalg.det(V)))) * np.exp(float((z-zeta).T * np.linalg.inv(V) * (z-zeta)))
    t2 = time.time()
    print(t2-t1)
    return sum

def lnprior(theta):
    alpha, beta, var = theta
    if -10.0 < alpha < 10.0 and -1.0 < beta < 1.0 and 0 < var < 10.0:
        return 0.0
    return -np.inf

def lnprob(theta, psi, Mvir, Mvir_err, cvir, cvir_err, redshift):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, psi, Mvir, Mvir_err, cvir, cvir_err, redshift)

nll = lambda *args: -lnlike(*args)



    
