import numpy as np
from scipy.misc import derivative
import MConvert_Personal as mc
import ipdb

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

def normalize_uncertainty(m,m_p,m_m,c,c_p,c_m,lamb=0.75):
    # making sure mass anc concentration are there
    assert m not in [u'TBD',u'nan'], "Mass is not present, or has yet to be determined..."
    assert c not in [u'TBD',u'nan'], "Concentration is not present, or has yet to be determined..."

    # making sure that m_p and c_p are always positive numbers
    if type(m_p) is float and np.isnan(m_p) is False:
        assert m_p >= 0, "Upper bound on mass is negative..."
    if type(c_p) is float and np.isnan(c_p) is False:
        assert c_p >= 0, "Upper bound on conc is negative..."
    # making sure that m_m and c_m are always negative numbers
    if type(m_m) is float and np.isnan(m_m) is False:
        assert m_m <= 0, "Lower bound on mass is positive..."
    if type(c_m) is float and np.isnan(c_m) is False:
        assert c_m <= 0, "Lower bound on conc is positive..."
    
    # determining if uncertainties are symmetric
    mass_bounds = False
    mass_is_sym = False
    conc_bounds = False
    conc_is_sym = False
    if m_p not in [u'TBD',u'nan',u'infty'] and m_m not in [u'TBD',u'nan']:
        mass_bounds = True
        if abs(m_p) == abs(m_m):
            mass_is_sym = True
    if c_p not in [u'TBD',u'nan',u'infty'] and c_m not in [u'TBD',u'nan']:
        conc_bounds = True
        if abs(c_p) == abs(c_m):
            conc_is_sym = True

    # symmetrization process according to D'Agostini, (using by default 75% of the
    # shift in the difference between the upper and lower bounds).
    if mass_bounds is True and mass_is_sym is False:
        m_sigtheta = (m_p + abs(m_m))/2.0
        m_expectedtheta = m + lamb*(m_p+m_m)
    if conc_bounds is True and conc_is_sym is False:
        c_sigtheta = (c_p + abs(c_m))/2.0
        c_expectedtheta = c + lamb*(c_p+c_m)

    # if mass uncertainties are there but conc. are not, adopt same symmetric
    # fractional uncertainty. 
    if mass_bounds is True and conc_bounds is False:
        if mass_is_sym is True:
            frac = abs(m_p)/m
            c_p = frac*c
            c_m = -1*frac*c
        elif mass_is_sym is False:
            frac = m_sigtheta/m_expectedtheta
            c_p = frac*c
            c_m = -1*frac*c
    # if mass uncertainties are there but conc. are not, adopt same symmetric
    # fractional uncertainty.
    if conc_bounds is True and mass_bounds is False:
        if conc_is_sym is True:
            frac = abs(c_p)/c
            m_p = frac*m
            m_m = -1*frac*m
        elif conc_is_sym is False:
            frac = c_sigtheta/c_expectedtheta
            m_p = frac*m
            m_m = -1*frac*m

    ## returning values ##

    # both uncertainties for both mass/conc
    # if both uncertainties are present for mass and conc, and they're both symmetric
    if mass_bounds is True and conc_bounds is True and mass_is_sym is True and conc_is_sym is True:
        return m,m_p,m_m,c,c_p,c_m
    # if both uncertainties are present for mass anc conc, and they're both asymmetric
    if mass_bounds is True and conc_bounds is True and mass_is_sym is False and conc_is_sym is False:
        return m_expectedtheta,m_sigtheta,-1*m_sigtheta,c_expectedtheta,c_sigtheta,-1*c_sigtheta
    # if both uncertainties are present for mass anc conc, and only mass is asymmetric
    if mass_bounds is True and conc_bounds is True and mass_is_sym is False and conc_is_sym is True:
        return m_expectedtheta,m_sigtheta,-1*m_sigtheta,c,c_p,c_m
    # if both uncertainties are present for mass anc conc, and only conc is asymmetric
    if mass_bounds is True and conc_bounds is True and mass_is_sym is True and conc_is_sym is False:
        return m,m_p,m_m,c_expectedtheta,c_sigtheta,-1*c_sigtheta

    # both uncertainties for just one of mass/conc
    # if mass uncertainties are present, and they're either symmetric or asymmetric
    if mass_bounds is True and conc_bounds is False:
        if mass_is_sym is True:
            return m,m_p,m_m,c,c_p,c_m
        elif mass_is_sym is False:
            return m_expectedtheta,m_sigtheta,-1*m_sigtheta,c,c_p,c_m
    # if conc uncertainties are present, and they're either symmetric or asymmetric
    if mass_bounds is False and conc_bounds is True:
        if conc_is_sym is True:
            return m,m_p,m_m,c,c_p,c_m
        elif conc_is_sym is False:
            return m,m_p,m_m,c_expectedtheta,c_sigtheta,-1*c_sigtheta

    # both uncertainties for mass/conc are missing; need to decide on an observationally motived way of assigning uncertainties
    if mass_bounds is False and conc_bounds is False:
        return m,0.5*m,-0.5*m,c,0.5*c,-0.5*c

def propagate_A_uncertainty(m,m_err,b,b_err,M_star=1.3e13/0.7):
    '''
    A function which propagates uncertainty through the relation for the normalization of the c-M model.
    '''
    A_best = 10**(b+m*np.log10(M_star))
    dAdb = np.log(10)*np.exp(np.log(10)*(b+m*np.log10(M_star)))
    dAdm = np.log10(M_star)*np.log(10)*np.exp(np.log(10)*(b+m*np.log10(M_star)))
    A_err = abs(dAdm)*m_err + abs(dAdb)*b_err
    return A_best,A_err

if __name__ == "__main__":

    '''
    # Converting measurements from 200 to virial
    z = 0.1114
    m200 = 7
    m200_err = 4.4
    c200 = 1.8
    c200_err = 1.5
    deltavir = mc.DeltaFinder(0.3,0.7,z)
    
    print "Mvir = " + "{}".format(mc.Mconvert(m200,200,deltavir,c200))
    print "cvir = " + "{}".format(mc.Cconvert(m200,200,deltavir,c200))
    print "Mvir_err = " + "{}".format(propagate_mass_uncertainty_independent(m200,m200_err,c200,c200_err,200,deltavir))
    print "cvir_err = " + "{}".format(propagate_conc_uncertainty(m200,m200_err,c200,c200_err,200,deltavir))
    '''
    ipdb.set_trace()
    #normalize(3.4,2.1,-2.1,7.4,1.9,-1.9)
