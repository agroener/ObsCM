import numpy as np
from scipy.misc import derivative



#################################
###          Functions        ###
#################################

def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = 1e-6)

def propagate_mass_uncertainty_independent(m_orig_best,m_orig_err,c_orig_best,c_orig_err,delta_orig,delta_new): # mass and concentration are part of the conversion process for mass
    partial_mass = partial_derivative(mc.Mconvert,var=0,point=[m_orig_best,delta_orig,delta_new,c_orig_best])
    partial_conc = partial_derivative(mc.Mconvert,var=3,point=[m_orig_best,delta_orig,delta_new,c_orig_best])
    uncertainty = np.sqrt((partial_mass*m_orig_err)**2 + (partial_conc*c_orig_err)**2)
    return uncertainty

def propagate_mass_uncertainty_notindependent(m_orig_best,m_orig_err,c_orig_best,c_orig_err,delta_orig,delta_new): # mass and concentration are part of the conversion process for mass
    partial_mass = partial_derivative(mc.Mconvert,var=0,point=[m_orig_best,delta_orig,delta_new,c_orig_best])
    partial_conc = partial_derivative(mc.Mconvert,var=3,point=[m_orig_best,delta_orig,delta_new,c_orig_best])
    uncertainty = abs(partial_mass*m_orig_err) + abs(partial_conc*c_orig_err)
    return uncertainty

def propagate_conc_uncertainty(m_orig_best,m_orig_err,c_orig_best,c_orig_err,delta_orig,delta_new): # mass is truly not part of the conversion process for concentration
    partial_conc = partial_derivative(mc.Cconvert,var=3,point=[m_orig_best,delta_orig,delta_new,c_orig_best])
    uncertainty = abs(partial_conc*c_orig_err)
    return uncertainty
