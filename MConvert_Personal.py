from __future__ import division # float division by default
import numpy as np
import ipdb

def Mconvert(Mold,DeltaOld,DeltaNew,ConcentOld):
    rat = DeltaOld/DeltaNew
    fval = np.log(1+ConcentOld)-ConcentOld/(1+ConcentOld)
    fval = fval/rat/(ConcentOld**3)
    Mconverted = Mold * ((ConcentOld*convinv(fval))**(-3)/rat)
    return Mconverted

def Cconvert(Mold,DeltaOld,DeltaNew,ConcentOld):
    rat = DeltaOld/DeltaNew
    fval = np.log(1+ConcentOld)-ConcentOld/(1+ConcentOld)
    fval = fval/rat/(ConcentOld**3)
    Cconverted = (convinv(fval))**(-1)
    return Cconverted

def convinv(f):
    a2 = 0.5116
    a3 = -1.285/3.
    a4 = -3.13e-3
    a5 = -3.52e-5
    p = a3+a4*np.log(f)+a5*(np.log(f))**2
    conv_inv = a2*f**(2*p)+(3./4.)**2
    conv_inv = 1./np.sqrt(conv_inv) + 2*f
    return conv_inv

def concentration_default(M,a):
    '''
    Default to Bullock et al in an LCDM sigma8=0.8 cosmology; assumes the mass, M, is in h^-1 M_sun
    '''
    Mstar = 2.77E12  # Mstar in h^-1 M_sun
    concent = 9.*a*(M/Mstar)**(-0.13)
    return concent

def DeltaFinder(Omega_m0,Omega_L0,z): # for flat universes only
    Omega_R0 = 1 - Omega_m0 - Omega_L0
    Ez2 = Omega_m0*(1+z)**3 + Omega_R0*(1+z)**2 + Omega_L0
    Omega = Omega_m0*(1+z)**3 / Ez2
    x = Omega - 1
    Delta = 18*np.pi**2 + 82*x - 39*x**2 
    return Delta

def DeltaFinder_Alt(Omega_m0,Omega_L0,z): # for flat universes only
    Omega_R0 = 1 - Omega_m0 - Omega_L0
    Ez2 = Omega_m0*(1+z)**3 + Omega_R0*(1+z)**2 + Omega_L0
    Omega = Omega_m0*(1+z)**3 / Ez2
    x = Omega - 1
    Delta = 18*np.pi**2 + 82*x - 32*x**2 # last term is supposed to be 39*x^2; trying out 32 for shits and giggles
    return Delta

def DeltaFinder2(Omega_m0,Omega_L0,z): # Hu & Kravtsov
    Omega_R0 = 1 - Omega_m0 - Omega_L0
    Ez2 = Omega_m0*(1+z)**3 + Omega_R0*(1+z)**2 + Omega_L0
    Omega = Omega_m0*(1+z)**3/Ez2
    x = Omega_m0 - 1
    Delta = (18*np.pi**2 + 82*x - 39*x**2) / (1+x)
    return Delta

def Dialog():
    
    #Setup
    print "Original Halo Mass?"
    Mold = raw_input('--> ')
    print "Original Concentration (ENTER for default)?"
    ConcentOld = raw_input('--> ')
    if not ConcentOld:
        ConcentOld = concentration_default(float(Mold),1.0)
        print "Repalcing with c={}".format(ConcentOld)
    print "Original Spherical Overdensity?"
    DeltaOld = raw_input('--> ')
    print 'Target Spherical Overdensity? (ENTER for over-density assistance)'
    DeltaNew = raw_input('--> ')
    if not DeltaNew:
        DeltaNew = Dialog2()
        print "DeltaNew = {}".format(DeltaNew)
 
    # Calculations and Output
    #pdb.set_trace()
    print "M_target = {0:.7}".format(Mconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld)))
    print "c_target = {0:.7}".format(Cconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld)))

    print "Repeat (yes/no)?"
    response = raw_input('--> ')
    if response == 'yes':
        return True
    elif response == 'no':
        return False

def Dialog2():
    print "   Omega_m_0?"
    Omega_m_0 = raw_input('   --> ')
    print "   Omega_L_0?"
    Omega_L_0 = raw_input('   --> ')
    print "   Redshift?"
    z = raw_input('   --> ')
    DeltaNew = DeltaFinder(float(Omega_m_0),float(Omega_L_0),float(z)) 
    return DeltaNew

def NoDialogToVir(Mold,ConcentOld,DeltaOld,redshift,Omega_m_0=0.3,Omega_L_0=0.7):
    DeltaNew = DeltaFinder(float(Omega_m_0),float(Omega_L_0),float(redshift))
    M_target = Mconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld))
    c_target = Cconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld))
    return M_target, c_target
    
def NoDialogFromVir(Mold,ConcentOld,DeltaNew,redshift,Omega_m_0=0.3,Omega_L_0=0.7):
    DeltaOld = DeltaFinder(float(Omega_m_0),float(Omega_L_0),float(redshift))
    M_target = Mconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld))
    c_target = Cconvert(float(Mold),float(DeltaOld),float(DeltaNew),float(ConcentOld))
    return M_target, c_target


if __name__ == '__main__':
    repeat = True
    while repeat:
        repeat = Dialog()

        
