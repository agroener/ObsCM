from astropy.io.ascii import read
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.misc import derivative
from scipy.optimize import leastsq
import numpy as np
import collections
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

# All Other Functions

def chi2(x,y,sigx,sigy,sig,m,b,N):
    sigtot2=sigy**2+sig**2
    sigx2=sigx**2
    x0=(m*y*sigx2-m*b*sigx2+x*sigtot2)/(m**2*sigx2+sigtot2)
    y0=m*x0+b
    chi2=sum((x-x0)**2/(sigx2+1.e-6)+(y-y0)**2/sigtot2)#/(N-3)
    return chi2

def partial_derivative(func, var=0, point=[]):
    args = point[:]
    def wraps(x):
        args[var] = x
        return func(*args)
    return derivative(wraps, point[var], dx = 1e-6)

def second_partial_derivative(func, var=0, point=[]):
    '''
    This function only works if taking second derivatives with respect
    to one variable. This does not work for arbitrary mixed partials.
    '''
    dx = 1e-6
    fprime = partial_derivative(func, var=var, point=point)
    point[var] += dx
    fprimep = partial_derivative(func, var=var, point=point)
    return (fprimep-fprime)/dx

def second_partial_derivative_mixed(func, var=[5,6], point=()):
    '''
    Somewhat specific to this problem, this function computes the
    second order mixed partial derivatives (between m, and b).
    '''
    delta = 1.e-1
    if var == [5,6]:
        # first calculate the various combinations of things
        x,y,sigx,sigy,sig,m,b,N = point
        fpp = func(x,y,sigx,sigy,sig,m+delta,b+delta,N)
        fmm = func(x,y,sigx,sigy,sig,m-delta,b-delta,N)
        fmp = func(x,y,sigx,sigy,sig,m-delta,b+delta,N)
        fpm = func(x,y,sigx,sigy,sig,m+delta,b-delta,N)
        return (1./(4*(delta**2)))*(fpp+fmm-fpm-fmp)
    else:
        print("Need to implement partial derivative...")


def startup(fname=None):
    if fname is not None:
        data=read(fname,delimiter=',')
        cl = data['col1'].data
        x = data['col2'].data
        y = data['col3'].data
        sigx = data['col4'].data
        sigy = data['col5'].data
        return x,y,sigx,sigy,cl
    else:
        print("Filename must be specified...")
        return

def startup_sims():
    try:
        raw_data = read('/Users/groenera/Desktop/Dropbox/Private/Research/GroupMeetings/Meeting#60/concs_m200s_data.txt',delimiter=',',guess=False)
    except:
        print("Cannot find datafile...")
        return

    raw_concs = raw_data['c200']
    raw_fof = raw_data['mfof']
    raw_masses = raw_data['m200']

    high_indices = [i for i in range(len(raw_fof)) if raw_fof[i] > 1.2e14]
    h_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in high_indices]
    h_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in high_indices]

    low_indices = [i for i in range(len(raw_masses)) if raw_fof[i] < 2.6e13]
    l_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in low_indices]
    l_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in low_indices]

    med_indices = [i for i in range(len(raw_masses)) if i not in high_indices and i not in low_indices]
    m_concs = [raw_concs[i] for i in range(len(raw_concs)) if i in med_indices]
    m_masses = [raw_masses[i] for i in range(len(raw_masses)) if i in med_indices]

    #plt.scatter(h_masses,h_concs,color='red')
    #plt.scatter(m_masses,m_concs,color='green')
    #plt.scatter(l_masses,l_concs,color='blue')
    #plt.xscale('log')
    #plt.show()

    return l_concs,l_masses,m_concs,m_masses,h_concs,h_masses

def startup_bootstrap(fname=None,handlerepeats=True,method=None):
    # first load in all of the data
    if handlerepeats is False:
        if fname is not None:
            data=read(fname,delimiter=',')
            cl = data['col1'].data
            x = data['col2'].data
            y = data['col3'].data
            sigx = data['col4'].data
            sigy = data['col5'].data
        else:
            print("Filename must be specified...")
            return
    elif handlerepeats is True:
        if method is not None:
            x,y,sigx,sigy,cl=discover_repeats(method=method)
        else:
            print("Method must be specified in startup_bootstrap if handlerepeats is True...")
    # now sample with replacement
    x_bs = np.random.choice(x,replace=True,size=len(x))
    bs_indices = [np.where(x == x_bs[i])[0][0] for i in range(len(x_bs))]
    y_bs = np.array([y[i] for i in bs_indices])
    sigx_bs = np.array([sigx[i] for i in bs_indices])
    sigy_bs = np.array([sigy[i] for i in bs_indices])
    cl_bs = np.array([cl[i] for i in bs_indices])
    return x_bs,y_bs,sigx_bs,sigy_bs,cl

def startup_bootstrap_all(handlerepeats=True):
    # first load in all of the data
    fname_list = ['CM_data.txt','LOSVD_data.txt','X-ray_data.txt',
                      'WL_data.txt','WL+SL_data.txt','SL_data.txt']
    if handlerepeats is False:
        cl,x,y,sigx,sigy = ([],[],[],[],[])
        if f in fname_list:
            data=read(f,delimiter=',')
            cl.append(data['col1'].data)
            x.append(data['col2'].data)
            y.append(data['col3'].data)
            sigx.append(data['col4'].data)
            sigy.append(data['col5'].data)
    elif handlerepeats is True:
        x,y,sigx,sigy,cl,methods=discover_repeats_all(fname_list=fname_list)
    # now sample with replacement
    x_bs = np.random.choice(x,replace=True,size=len(x))
    bs_indices = [np.where(x == x_bs[i])[0][0] for i in range(len(x_bs))]
    y_bs = np.array([y[i] for i in bs_indices])
    sigx_bs = np.array([sigx[i] for i in bs_indices])
    sigy_bs = np.array([sigy[i] for i in bs_indices])
    cl_bs = np.array([cl[i] for i in bs_indices])
    return x_bs,y_bs,sigx_bs,sigy_bs,cl

def steepest_decent(x,y,sigx,sigy,sig,m,b,N,alpha=0.8,tol=1e-6):

    converged = False

    while converged is False:

        # Calculate the chi2 of the current point
        chi0 = chi2(x,y,sigx,sigy,sig,m,b,N)
        # Calculate first derivatives at point
        dchi2_db = partial_derivative(chi2,var=6,point=[x,y,sigx,sigy,sig,m,b,N])
        dchi2_dm = partial_derivative(chi2,var=5,point=[x,y,sigx,sigy,sig,m,b,N])
        # Calculate second derivatives at point
        d2chi2_db2 = second_partial_derivative(chi2,var=6,point=[x,y,sigx,sigy,sig,m,b,N])
        d2chi2_dm2 = second_partial_derivative(chi2,var=5,point=[x,y,sigx,sigy,sig,m,b,N])
        delta_b=-alpha*dchi2_db/d2chi2_db2
        delta_m=-alpha*dchi2_dm/d2chi2_dm2
        # Check to see if steps in either direction are better
        chi_new_b = chi2(x,y,sigx,sigy,sig,m,b+delta_b,N)
        chi_new_m = chi2(x,y,sigx,sigy,sig,m+delta_m,b,N)
        # Check for convergence; else take that step
        if abs(chi_new_b-chi0) <= tol and abs(chi_new_m-chi0):
            print("Chi2_red when converged: {}".format(chi0/(N-2)))# N-2 since only fitting for m and b currently
            converged = True
        else:
            b = b + delta_b
            m = m + delta_m

    return m,b

def fit(method='X-ray', plot=False, savefig=False):
    # Load the data from file
    if method in ['X-ray','x-ray','xray']:
        pl_col = 'green'
        method_title = 'X-ray'
        filename='{}_data.txt'.format('X-ray')
    elif method in ['WL','wl']:
        pl_col = 'purple'
        method_title = 'WL'
        filename='{}_data.txt'.format('WL')
    elif method in ['SL','sl']:
        pl_col = 'red'
        method_title = 'SL'
        filename='{}_data.txt'.format('SL')
    elif method in ['WL+SL','wl+sl']:
        pl_col = 'black'
        method_title = 'WL+SL'
        filename='{}_data.txt'.format('WL+SL')
    elif method in ['CM','cm']:
        pl_col = 'blue'
        method_title = 'CM'
        filename='{}_data.txt'.format('CM')
    elif method in ['LOSVD','losvd']:
        pl_col = 'orange'
        method_title = 'LOSVD'
        filename='{}_data.txt'.format('LOSVD')
    else:
        print("Method is undefined...")
        return

    intro = "FITTING METHOD: {}".format(method_title)
    print('*'*len(intro))
    print(intro)
    print('*'*len(intro))

    x_old,y_old,sigx_old,sigy_old,cl_old = startup(fname=filename) # doesn't take repeat measurements into account; but need it for original sample size
    x,y,sigx,sigy,cl=discover_repeats(method=method)
    x,y,sigx,sigy = (np.array(x),np.array(y),np.array(sigx),np.array(sigy))
    xmax,xmin = max(x),min(x)

    print("Number of measurements used: {}".format(len(x_old)))
    print("Number of unique clusters: {}".format(len(x)))

    # Fitting as if there are no uncertainties
    N=len(x)
    Sxy=sum(x*y)
    Sx=sum(x)
    Sy=sum(y)
    Sxx=sum(x*x)
    m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
    b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
    sig=np.sqrt(np.std(y-(m1*x+b1))**2-np.mean(sigy**2))
    print("Linear model (assuming no uncertainties): m: {}, b: {}, sig: {} ".format(m1,b1,sig))

    # Using chi-squared fitting routine
    ## (not completely the correct thing to do,
    ##  but a back of the envelope method)
    m2,b2 = steepest_decent(x,y,sigx,sigy,sig,m1,b1,N,alpha=0.75,tol=1.e-6)

    # Calculating second order partial derivatives for error of the fit
    F11 = second_partial_derivative(chi2,var=5,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F22 = second_partial_derivative(chi2,var=6,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F12 = second_partial_derivative_mixed(chi2,var=[5,6],point=(x,y,sigx,sigy,sig,m2,b2,N))
    F21 = F12
    # Error in fit
    sigm=1./np.sqrt(F11)
    sigb=1./np.sqrt(F22)

    print("Linear model (with uncertainties): m: {} +/- {}, b: {} +/- {}".format(m2,sigm,b2,sigb))


    # Plotting the best fit line over the data
    if plot is True:
        plt.figure(figsize=(8,8))
        plt.title(method_title)
        plt.errorbar(x,y,xerr=sigx,yerr=sigy,fmt='o',color=pl_col)
        x0 = np.array([13.0,17.5])
        y0=m2*x0+b2
        plt.plot(x0,y0,linewidth=3,color=pl_col)
        plt.fill_between(x0,[(m2+sigm)*i+(b2+sigb+sig) for i in x0],[(m2-sigm)*i+(b2-sigb-sig) for i in x0],alpha=0.25,color=pl_col)
        plt.xlabel(r'$\mathrm{\log{ M_{vir}/M_{\odot}}}$',fontsize=18)
        plt.ylabel(r'$\mathrm{\log{ \, c_{vir} \, (1+z) }}$',fontsize=18)
        plt.xlim(13.0,17.5)
        plt.ylim(-1.0,2.5)
        if savefig is True:
            plt.savefig('{}_linearmodel_witherror.png'.format(filename.split('_')[0]))

    # Calculating and plotting error ellipse
    detF=F11*F22-F12*F12
    Qxx=F22/detF
    Qyy=F11/detF
    Qxy=-F12/detF

    theta=np.arctan(2*Qxy/(Qxx-Qyy+1e-9))/2.
    a1=np.sqrt(Qxx*pow(np.cos(theta),2)+Qyy*pow(np.sin(theta),2)+2*Qxy*np.sin(theta)*np.cos(theta))
    b1=np.sqrt(Qxx*pow(np.sin(theta),2)+Qyy*pow(np.cos(theta),2)-2*Qxy*np.sin(theta)*np.cos(theta))

    if plot is True:
        fig = plt.figure(figsize=(8,8))
        plt.title(method_title)
        ax = fig.add_subplot(111)
        ell=Ellipse([m2,b2],width=2*a1,height=2*b1,angle=theta*57.3)
        ell2=Ellipse([m2,b2],width=a1,height=b1,angle=theta*57.3)
        ax.add_artist(ell)
        ell.set_facecolor('g')
        ax.add_artist(ell2)
        ell2.set_facecolor('y')
        ax.plot(m2,b2,'+')
        plt.xlabel('m')
        plt.ylabel('b')
        ax.set_xlim(m2-30*a1,m2+30*a1)
        ax.set_ylim(b2-300*a1,b2+300*a1)
        if savefig is True:
            plt.savefig('{}_errorellipse.png'.format(filename.split('_')[0]))
        plt.show()

    return m2,sigm,b2,sigb,sig,(a1,b1,theta),(xmax,xmin)

def fit_bootstrap(method=None, witherrors=True, nsamples=100):
    # Load the data from file
    if method in ['X-ray','x-ray','xray']:
        pl_col = 'green'
        method_title = 'X-ray'
        filename='{}_data.txt'.format('X-ray')
    elif method in ['WL','wl']:
        pl_col = 'purple'
        method_title = 'WL'
        filename='{}_data.txt'.format('WL')
    elif method in ['SL','sl']:
        pl_col = 'red'
        method_title = 'SL'
        filename='{}_data.txt'.format('SL')
    elif method in ['WL+SL','wl+sl']:
        pl_col = 'black'
        method_title = 'WL+SL'
        filename='{}_data.txt'.format('WL+SL')
    elif method in ['CM','cm']:
        pl_col = 'blue'
        method_title = 'CM'
        filename='{}_data.txt'.format('CM')
    elif method in ['LOSVD','losvd']:
        pl_col = 'orange'
        method_title = 'LOSVD'
        filename='{}_data.txt'.format('LOSVD')
    else:
        print("Method is undefined...")
        return

    if witherrors is False:
        m1_list, b1_list, sig_list = ([],[],[])
        for i in range(nsamples):
            x_bs,y_bs,sigx_bs,sigy_bs,cl = startup_bootstrap(fname=filename,method=method)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=np.sqrt(abs(np.std(y_bs-(m1*x_bs+b1))**2-np.mean(sigy_bs**2)))
            m1_list.append(m1)
            b1_list.append(b1)
            sig_list.append(sig)
        return m1_list, b1_list, sig_list, sig

    if witherrors is True:
        m2_list, b2_list = ([],[])
        for i in range(nsamples):
            x_bs,y_bs,sigx_bs,sigy_bs,cl = startup_bootstrap(fname=filename,method=method)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=np.sqrt(abs(np.std(y_bs-(m1*x_bs+b1))**2-np.mean(sigy_bs**2))) #abs needed for cases where second term is larger
            m2,b2 = steepest_decent(x_bs,y_bs,sigx_bs,sigy_bs,sig,m1,b1,N,alpha=0.75,tol=1.e-6)
            m2_list.append(m2)
            b2_list.append(b2)
        return m2_list, b2_list, sig
    

## Consider doing a bootstrap estimate of intrinsic scatter
def fit_bootstrap_allmethods(witherrors=True, nsamples=100):
    if witherrors is False:
        m1_list, b1_list, sig_list = ([],[],[])
        for i in range(nsamples):
            x_bs,y_bs,sigx_bs,sigy_bs,cl = startup_bootstrap_all(handlerepeats=True)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=np.sqrt(abs(np.std(y_bs-(m1*x_bs+b1))**2-np.mean(sigy_bs**2)))
            m1_list.append(m1)
            b1_list.append(b1)
            sig_list.append(sig)
        return m1_list, b1_list, sig_list, sig

    if witherrors is True:
        m2_list, b2_list = ([],[])
        for i in range(nsamples):
            x_bs,y_bs,sigx_bs,sigy_bs,cl = startup_bootstrap_all(handlerepeats=True)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=np.sqrt(abs(np.std(y_bs-(m1*x_bs+b1))**2-np.mean(sigy_bs**2))) #abs needed for cases where second term is larger
            m2,b2 = steepest_decent(x_bs,y_bs,sigx_bs,sigy_bs,sig,m1,b1,N,alpha=0.75,tol=1.e-6)
            m2_list.append(m2)
            b2_list.append(b2)
        return m2_list, b2_list, sig
    
def fit_sims(los_angle=None,scale=None):
    l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()
    # Fitting with no uncertainties
    if los_angle is None:
        l_concs = np.log10(np.array(l_concs))
        l_masses = np.log10(np.array(l_masses))
        m_concs = np.log10(np.array(m_concs))
        m_masses = np.log10(np.array(m_masses))
        h_concs = np.log10(np.array(h_concs))
        h_masses = np.log10(np.array(h_masses))
        x = np.array(list(l_masses)+list(m_masses)+list(h_masses))
        y = np.array(list(l_concs)+list(m_concs)+list(h_concs))
    elif los_angle in ['Aligned','aligned']:
        if scale in ['r200','R200','full']:
            low_shapes = [np.random.normal(loc=0.64,scale=0.12) for i in range(len(l_concs))]
            med_shapes = [np.random.normal(loc=0.59,scale=0.115) for i in range(len(m_concs))]
            high_shapes = [np.random.normal(loc=0.48,scale=0.115) for i in range(len(h_concs))]
        elif scale in ['half','0.5*r200','0.5*R200']:
            low_shapes = [np.random.normal(loc=0.60,scale=0.125) for i in range(len(l_concs))]
            med_shapes = [np.random.normal(loc=0.555,scale=0.125) for i in range(len(m_concs))]
            high_shapes = [np.random.normal(loc=0.49,scale=0.095) for i in range(len(h_concs))]
        else:
            print("Invalid input for variable scale...")
            return
        # now recompute concentration as if aligned along the l.o.s.
        y_low = [conc_finder_pro(l_concs[i],[0],low_shapes[i])[0][0] for i in range(len(l_concs))]
        y_med = [conc_finder_pro(m_concs[i],[0],med_shapes[i])[0][0] for i in range(len(m_concs))]
        y_high = [conc_finder_pro(h_concs[i],[0],high_shapes[i])[0][0] for i in range(len(h_concs))]
        y = np.log10(np.array(y_low+y_med+y_high))
        x = np.log10(np.array(list(l_masses)+list(m_masses)+list(h_masses)))

        #f, axarr = plt.subplots(2, sharex=True)
        #axarr[0].hist(l_concs,bins=30,color='blue',histtype='stepfilled',alpha=0.5,normed=True)
        #axarr[0].hist(m_concs,bins=25,color='green',histtype='stepfilled',alpha=0.5,normed=True)
        #axarr[0].hist(m_concs,bins=5,color='red',histtype='stepfilled',alpha=0.5,normed=True)
        #axarr[1].hist(y_low,bins=30,color='blue',histtype='stepfilled',alpha=0.5,normed=True)
        #axarr[1].hist(y_med,bins=25,color='green',histtype='stepfilled',alpha=0.5,normed=True)
        #axarr[1].hist(y_high,bins=5,color='red',histtype='stepfilled',alpha=0.5,normed=True)
        #plt.show()
        #ipdb.set_trace()

    N=len(x)
    Sxy=sum(x*y)
    Sx=sum(x)
    Sy=sum(y)
    Sxx=sum(x*x)
    m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
    b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
    sig=np.std(y-(m1*x+b1))
    print("Linear model (assuming no uncertainties): m: {}, b: {}, sig: {} ".format(m1,b1,sig))
    return m1,b1,sig

def fit_all(plot=False, savefig=False, plotwitherrorbars = False):
    # Load the data from file
    pl_col_xray = 'green'
    filename_xray='{}_data.txt'.format('X-ray')
    pl_col_wl = 'purple'
    filename_wl='{}_data.txt'.format('WL')
    pl_col_sl = 'red'
    filename_sl='{}_data.txt'.format('SL')
    pl_col_wlsl = 'black'
    filename_wlsl='{}_data.txt'.format('WL+SL')
    pl_col_cm = 'blue'
    filename_cm='{}_data.txt'.format('CM')
    pl_col_losvd = 'orange'
    filename_losvd='{}_data.txt'.format('LOSVD')

    intro = "FITTING METHOD: FULL OBSERVATIONAL SAMPLE"
    print('*'*len(intro))
    print(intro)
    print('*'*len(intro))

    # doesn't take repeat measurements into account; but need it to print the original sample size to the terminal
    tot = 0
    for i in [filename_xray,filename_wl,filename_sl,filename_wlsl,filename_cm,filename_losvd]:
        x_old,y_old,sigx_old,sigy_old,cl_old = startup(fname=i)
        tot = tot + len(x_old)
    
    # new function which discovers repeats of clusters across any and all methods
    x,y,sigx,sigy,uniques,methods = discover_repeats_all(
        fname_list=[filename_xray,filename_wl,filename_sl,
                    filename_wlsl,filename_cm,filename_losvd])
    x,y,sigx,sigy = (np.array(x),np.array(y),np.array(sigx),np.array(sigy))
    xmax,xmin = max(x),min(x)
    
    print("Number of measurements used: {}".format(tot))
    print("Number of unique clusters: {}".format(len(x)))

    # Fitting as if there are no uncertainties
    N=len(x)
    Sxy=sum(x*y)
    Sx=sum(x)
    Sy=sum(y)
    Sxx=sum(x*x)
    m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
    b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
    sig=np.sqrt(np.std(y-(m1*x+b1))**2-np.mean(sigy**2))
    print("Linear model (assuming no uncertainties): m: {}, b: {}, sig: {} ".format(m1,b1,sig))

    # Using chi-squared fitting routine
    ## (not completely the correct thing to do,
    ##  but a back of the envelope method)
    m2,b2 = steepest_decent(x,y,sigx,sigy,sig,m1,b1,N,alpha=0.75,tol=1.e-6)

    # Calculating second order partial derivatives for error of the fit
    F11 = second_partial_derivative(chi2,var=5,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F22 = second_partial_derivative(chi2,var=6,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F12 = second_partial_derivative_mixed(chi2,var=[5,6],point=(x,y,sigx,sigy,sig,m2,b2,N))
    F21 = F12
    # Error in fit
    sigm=1./np.sqrt(F11)
    sigb=1./np.sqrt(F22)

    print("Linear model (with uncertainties): m: {} +/- {}, b: {} +/- {}".format(m2,sigm,b2,sigb))
    
    # Plotting the best fit line over the data
    if plot is True:
        color_list = []
        plt.figure(figsize=(8,8),dpi=120)
        plt.title("All Methods")
        for i in range(len(x)):
            if methods[i] == 'X-ray':
                tmp_col = 'green'
            elif methods[i] == 'WL':
                tmp_col = 'purple'
            elif methods[i] == 'SL':
                tmp_col = 'red'
            elif methods[i] == 'WL+SL':
                tmp_col = 'black'
            elif methods[i] == 'CM':
                tmp_col = 'blue'
            elif methods[i] == 'LOSVD':
                tmp_col = 'orange'
            elif methods[i] == 'COMB':
                tmp_col = 'cyan'
            else:
                print("Undefined method found...")
                return
            if plotwitherrorbars is True:
                plt.errorbar(x[i],y[i],xerr=sigx[i],yerr=sigy[i],fmt='o',color=tmp_col)
            elif plotwitherrorbars is False:
                plt.scatter(x[i],y[i],color=tmp_col)
        x0 = np.array([13.0,17.5])
        y0=m2*x0+b2
        plt.plot(x0,y0,linewidth=3,color='yellow')
        plt.fill_between(x0,[(m2+sigm)*i+(b2+sigb+sig) for i in x0],[(m2-sigm)*i+(b2-sigb-sig) for i in x0],alpha=0.5,color='yellow')
        l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()
        plt.errorbar(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),yerr=np.log10(np.std(l_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,label='Groener and Goldberg (2014)',zorder=20)
        plt.scatter(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.errorbar(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),yerr=np.log10(np.std(m_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20)
        plt.scatter(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.errorbar(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),yerr=np.log10(np.std(h_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20)
        plt.scatter(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.xlabel(r'$\mathrm{\log{ M_{vir}/M_{\odot}}}$',fontsize=25)
        plt.ylabel(r'$\mathrm{\log{ \, c_{vir} \, (1+z) }}$',fontsize=25)
        plt.xlim(13.0,17.0)
        plt.ylim(-0.5,2.0)
        plt.legend(loc=0,numpoints=1,frameon=False,fontsize=12)
        if savefig is True:
            plt.savefig('AllMethodsWithSims_linearmodel_witherror.png')

    # Calculating and plotting error ellipse
    detF=F11*F22-F12*F12
    Qxx=F22/detF
    Qyy=F11/detF
    Qxy=-F12/detF

    theta=np.arctan(2*Qxy/(Qxx-Qyy+1e-9))/2.
    a1=np.sqrt(Qxx*pow(np.cos(theta),2)+Qyy*pow(np.sin(theta),2)+2*Qxy*np.sin(theta)*np.cos(theta))
    b1=np.sqrt(Qxx*pow(np.sin(theta),2)+Qyy*pow(np.cos(theta),2)-2*Qxy*np.sin(theta)*np.cos(theta))

    if plot is True:
        fig = plt.figure(figsize=(8,8))
        plt.title("All Methods")
        ax = fig.add_subplot(111)
        ell=Ellipse([m2,b2],width=2*a1,height=2*b1,angle=theta*57.3)
        ell2=Ellipse([m2,b2],width=a1,height=b1,angle=theta*57.3)
        ax.add_artist(ell)
        ell.set_facecolor('g')
        ax.add_artist(ell2)
        ell2.set_facecolor('y')
        ax.plot(m2,b2,'+')
        plt.xlabel('m')
        plt.ylabel('b')
        ax.set_xlim(m2-30*a1,m2+30*a1)
        ax.set_ylim(b2-300*a1,b2+300*a1)
        if savefig is True:
            plt.savefig('AllMethods_errorellipse.png')
        plt.show()

    return m2,sigm,b2,sigb,sig,(a1,b1,theta),(xmax,xmin)

def fit_multimeas_clusters(metho1=None, method2=None):
    # Load the data from file
    pl_col_xray = 'green'
    filename_xray='{}_data.txt'.format('X-ray')
    pl_col_wl = 'purple'
    filename_wl='{}_data.txt'.format('WL')
    pl_col_sl = 'red'
    filename_sl='{}_data.txt'.format('SL')
    pl_col_wlsl = 'black'
    filename_wlsl='{}_data.txt'.format('WL+SL')
    pl_col_cm = 'blue'
    filename_cm='{}_data.txt'.format('CM')
    pl_col_losvd = 'orange'
    filename_losvd='{}_data.txt'.format('LOSVD')

    intro = "FITTING METHOD: FULL OBSERVATIONAL SAMPLE"
    print('*'*len(intro))
    print(intro)
    print('*'*len(intro))

    # doesn't take repeat measurements into account; but need it for print the original sample size to the terminal
    tot = 0
    for i in [filename_xray,filename_wl,filename_sl,filename_wlsl,filename_cm,filename_losvd]:
        x_old,y_old,sigx_old,sigy_old,cl_old = startup(fname=i)
        tot = tot + len(x_old)

    # new function which discovers repeats of clusters across any and all methods
    x,y,sigx,sigy,uniques,methods = discover_repeats_all(
        fname_list=[filename_xray,filename_wl,filename_sl,
                    filename_wlsl,filename_cm,filename_losvd])
    x,y,sigx,sigy = (np.array(x),np.array(y),np.array(sigx),np.array(sigy))
    xmax,xmin = max(x),min(x)

    print("Number of measurements used: {}".format(tot))
    print("Number of unique clusters: {}".format(len(x)))

    # Fitting as if there are no uncertainties
    N=len(x)
    Sxy=sum(x*y)
    Sx=sum(x)
    Sy=sum(y)
    Sxx=sum(x*x)
    m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
    b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
    sig=np.sqrt(np.std(y-(m1*x+b1))**2-np.mean(sigy**2))
    print("Linear model (assuming no uncertainties): m: {}, b: {}, sig: {} ".format(m1,b1,sig))

    # Using chi-squared fitting routine
    ## (not completely the correct thing to do,
    ##  but a back of the envelope method)
    m2,b2 = steepest_decent(x,y,sigx,sigy,sig,m1,b1,N,alpha=0.75,tol=1.e-6)

    # Calculating second order partial derivatives for error of the fit
    F11 = second_partial_derivative(chi2,var=5,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F22 = second_partial_derivative(chi2,var=6,point=[x,y,sigx,sigy,sig,m2,b2,N])
    F12 = second_partial_derivative_mixed(chi2,var=[5,6],point=(x,y,sigx,sigy,sig,m2,b2,N))
    F21 = F12
    # Error in fit
    sigm=1./np.sqrt(F11)
    sigb=1./np.sqrt(F22)

    print("Linear model (with uncertainties): m: {} +/- {}, b: {} +/- {}".format(m2,sigm,b2,sigb))

    # Plotting the best fit line over the data
    if plot is True:
        color_list = []
        plt.figure(figsize=(8,8))
        plt.title("All Methods")
        for i in range(len(x)):
            if methods[i] == 'X-ray':
                tmp_col = 'green'
            elif methods[i] == 'WL':
                tmp_col = 'purple'
            elif methods[i] == 'SL':
                tmp_col = 'red'
            elif methods[i] == 'WL+SL':
                tmp_col = 'black'
            elif methods[i] == 'CM':
                tmp_col = 'blue'
            elif methods[i] == 'LOSVD':
                tmp_col = 'orange'
            elif methods[i] == 'COMB':
                tmp_col = 'cyan'
            else:
                print("Undefined method found...")
                return
            if plotwitherrorbars is True:
                plt.errorbar(x[i],y[i],xerr=sigx[i],yerr=sigy[i],fmt='o',color=tmp_col)
            elif plotwitherrorbars is False:
                plt.scatter(x[i],y[i],color=tmp_col)
        x0 = np.array([13.0,17.5])
        y0=m2*x0+b2
        plt.plot(x0,y0,linewidth=3,color='yellow')
        plt.fill_between(x0,[(m2+sigm)*i+(b2+sigb+sig) for i in x0],[(m2-sigm)*i+(b2-sigb-sig) for i in x0],alpha=0.5,color='yellow')
        l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()
        plt.errorbar(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),yerr=np.log10(np.std(l_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,label='Groener and Goldberg (2014)',zorder=20)
        plt.scatter(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.errorbar(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),yerr=np.log10(np.std(m_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20)
        plt.scatter(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.errorbar(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),yerr=np.log10(np.std(h_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20)
        plt.scatter(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k')
        plt.xlabel(r'$\mathrm{\log{ M_{vir}/M_{\odot}}}$',fontsize=18)
        plt.ylabel(r'$\mathrm{\log{ \, c_{vir} \, (1+z) }}$',fontsize=18)
        plt.xlim(13.0,17.0)
        plt.ylim(-0.5,2.0)
        plt.legend(loc=0,numpoints=1,frameon=False,fontsize=11)
        if savefig is True:
            plt.savefig('AllMethodsWithSims_linearmodel_witherror.png')

    # Calculating and plotting error ellipse
    detF=F11*F22-F12*F12
    Qxx=F22/detF
    Qyy=F11/detF
    Qxy=-F12/detF

    theta=np.arctan(2*Qxy/(Qxx-Qyy+1e-9))/2.
    a1=np.sqrt(Qxx*pow(np.cos(theta),2)+Qyy*pow(np.sin(theta),2)+2*Qxy*np.sin(theta)*np.cos(theta))
    b1=np.sqrt(Qxx*pow(np.sin(theta),2)+Qyy*pow(np.cos(theta),2)-2*Qxy*np.sin(theta)*np.cos(theta))

    if plot is True:
        fig = plt.figure(figsize=(8,8))
        plt.title("All Methods")
        ax = fig.add_subplot(111)
        ell=Ellipse([m2,b2],width=2*a1,height=2*b1,angle=theta*57.3)
        ell2=Ellipse([m2,b2],width=a1,height=b1,angle=theta*57.3)
        ax.add_artist(ell)
        ell.set_facecolor('g')
        ax.add_artist(ell2)
        ell2.set_facecolor('y')
        ax.plot(m2,b2,'+')
        plt.xlabel('m')
        plt.ylabel('b')
        ax.set_xlim(m2-30*a1,m2+30*a1)
        ax.set_ylim(b2-300*a1,b2+300*a1)
        if savefig is True:
            plt.savefig('AllMethods_errorellipse.png')
        plt.show()

    return m2,sigm,b2,sigb,sig,(a1,b1,theta),(xmax,xmin)


def do_bootstrap(method = 'X-ray', nsamples = 1000, witherrors = False, showplots=False, savepath=None):
    if method in ['X-ray','x-ray','xray']:
        pl_col = 'green'
        method_title = 'X-ray'
        filename='{}_data.txt'.format('X-ray')
    elif method in ['WL','wl']:
        pl_col = 'purple'
        method_title = 'WL'
        filename='{}_data.txt'.format('WL')
    elif method in ['SL','sl']:
        pl_col = 'red'
        method_title = 'SL'
        filename='{}_data.txt'.format('SL')
    elif method in ['WL+SL','wl+sl']:
        pl_col = 'black'
        method_title = 'WL+SL'
        filename='{}_data.txt'.format('WL+SL')
    elif method in ['CM','cm']:
        pl_col = 'blue'
        method_title = 'CM'
        filename='{}_data.txt'.format('CM')
    elif method in ['LOSVD','losvd']:
        pl_col = 'orange'
        method_title = 'LOSVD'
        filename='{}_data.txt'.format('LOSVD')
    else:
        print("Method is undefined...")
        return

    if savepath is None:
        path = ''
    else:
        path = savepath

    if witherrors is False:
        err_str = 'withouterr'
        m, b, sig = fit_bootstrap(method=method, witherrors=witherrors, nsamples=nsamples)
        m = [m[i] for i in range(len(m)) if not np.isnan(m[i])]
        b = [b[i] for i in range(len(b)) if not np.isnan(b[i])]
        sig = [sig[i] for i in range(len(sig)) if not np.isnan(sig[i])]
        text_height1 = 35
        text_height2 = 25

    elif witherrors is True:
        err_str = 'witherr'
        m, b, sig = fit_bootstrap(method=method, witherrors=witherrors, nsamples=nsamples)
        m = [m[i] for i in range(len(m)) if not np.isnan(m[i])]
        b = [b[i] for i in range(len(b)) if not np.isnan(b[i])]
        text_height1 = 3.5
        text_height2 = 2.5


    plt.figure()
    plt.hist(m,bins=30,color=pl_col,histtype='stepfilled',alpha=0.5)
    plt.title('Number of Bootstrap Samples: {}'.format(nsamples))
    plt.xlabel(r'$\mathrm{m}$',fontsize=18)
    plt.axvline(x=np.average(m),color='black',linewidth=3)
    plt.axvline(x=np.average(m)+np.std(m),linestyle='--',color='black',linewidth=3)
    plt.axvline(x=np.average(m)-np.std(m),linestyle='--',color='black',linewidth=3)
    plt.text(np.average(m)-3*np.std(m),text_height1,"mean: {:.4f}".format(np.average(m)))
    plt.text(np.average(m)-3*np.std(m),text_height2,"std.: {:.4f}".format(np.std(m)))
    plt.savefig(path+"{}_m_bootstrap_{}.png".format(method_title,err_str))

    plt.figure()
    plt.hist(b,bins=30,color=pl_col,histtype='stepfilled',alpha=0.5)
    plt.title('Number of Bootstrap Samples: {}'.format(nsamples))
    plt.xlabel(r'$\mathrm{b}$',fontsize=18)
    plt.axvline(x=np.average(b),color='black',linewidth=3)
    plt.axvline(x=np.average(b)+np.std(b),linestyle='--',color='black',linewidth=3)
    plt.axvline(x=np.average(b)-np.std(b),linestyle='--',color='black',linewidth=3)
    plt.text(np.average(b)-3*np.std(b),text_height1,"mean: {:.4f}".format(np.average(b)))
    plt.text(np.average(b)-3*np.std(b),text_height2,"std.: {:.4f}".format(np.std(b)))
    plt.savefig(path+"{}_b_bootstrap_{}.png".format(method_title,err_str))

    if witherrors is False:
        plt.figure()
        plt.hist(sig,bins=30,color=pl_col,histtype='stepfilled',alpha=0.5)
        plt.title('Number of Bootstrap Samples: {}'.format(nsamples))
        plt.xlabel(r'$\mathrm{\sigma}$',fontsize=18)
        plt.axvline(x=np.average(sig),color='black',linewidth=3)
        plt.axvline(x=np.average(sig)+np.std(sig),linestyle='--',color='black',linewidth=3)
        plt.axvline(x=np.average(sig)-np.std(sig),linestyle='--',color='black',linewidth=3)
        plt.text(np.average(sig)-3*np.std(sig),text_height1,"mean: {:.4f}".format(np.average(sig)))
        plt.text(np.average(sig)-3*np.std(sig),text_height2,"std.: {:.4f}".format(np.std(sig)))
        plt.savefig(path+"{}_sig_bootstrap_{}.png".format(method_title,err_str))

    if showplots is True:
        plt.show()

def discover_repeats(method=None, plotrepeats=False, plotdiff=False):

    if method in ['X-ray','x-ray','xray']:
        pl_col = 'green'
        method_title = 'X-ray'
        filename='{}_data.txt'.format('X-ray')
    elif method in ['WL','wl']:
        pl_col = 'purple'
        method_title = 'WL'
        filename='{}_data.txt'.format('WL')
    elif method in ['SL','sl']:
        pl_col = 'red'
        method_title = 'SL'
        filename='{}_data.txt'.format('SL')
    elif method in ['WL+SL','wl+sl']:
        pl_col = 'black'
        method_title = 'WL+SL'
        filename='{}_data.txt'.format('WL+SL')
    elif method in ['CM','cm']:
        pl_col = 'blue'
        method_title = 'CM'
        filename='{}_data.txt'.format('CM')
    elif method in ['LOSVD','losvd']:
        pl_col = 'orange'
        method_title = 'LOSVD'
        filename='{}_data.txt'.format('LOSVD')
    else:
        print("Method is undefined...")
        return

    x,y,sigx,sigy,clusters = startup(fname=filename)
    clusters = clusters.tolist()
    uniques = [i for i in set(clusters)]
    counts = [clusters.count(i) for i in uniques]
    xout = []
    yout = []
    sigxout = []
    sigyout = []
    for i in range(len(uniques)):
        if counts[i] == 1:
            tmp_index = clusters.index(uniques[i])
            xout.append(x[tmp_index])
            yout.append(y[tmp_index])
            sigxout.append(sigx[tmp_index])
            sigyout.append(sigy[tmp_index])
        else:
            tmp_indices = [j for j in range(len(clusters)) if clusters[j] == uniques[i]]
            xrep = [x[i] for i in tmp_indices]
            sigxrep = [sigx[i] for i in tmp_indices]
            xweights = [(1.0/sigxrep[i]**2) for i in range(len(xrep))]
            xnew = sum([xrep[i]*xweights[i] for i in range(len(xrep))])/sum(xweights)
            sigxnew = 1.0/np.sqrt(sum(xweights))
            xout.append(xnew)
            sigxout.append(sigxnew)
            yrep = [y[i] for i in tmp_indices]
            sigyrep = [sigy[i] for i in tmp_indices]
            yweights = [(1.0/sigyrep[i]**2) for i in range(len(yrep))]
            ynew = sum([yrep[i]*yweights[i] for i in range(len(yrep))])/sum(yweights)
            sigynew = 1.0/np.sqrt(sum(yweights))
            yout.append(ynew)
            sigyout.append(sigynew)
            if plotrepeats is True:
                plt.title("{}".format(uniques[i]))
                plt.xlabel(r"$\mathrm{\log M_{vir}/M_{\odot}}$",fontsize=18)
                plt.ylabel(r"$\mathrm{\log c_{vir} (1+z)}$",fontsize=18)
                plt.errorbar(xnew,ynew,yerr=[sigynew],xerr=[sigxnew],color='k')
                for i in range(len(tmp_indices)):
                    plt.errorbar(xrep[i],yrep[i],yerr=[sigyrep[i]],xerr=[sigxrep[i]],color=pl_col)
                plt.show()
                ipdb.set_trace()
    if plotdiff is True:
        f, axarr = plt.subplots(2, sharex=True)
        f.text(0.5, 0.01, r"$\mathrm{\log \, M_{vir}/M_{\odot}}$", fontsize=18, ha='center', va='center')
        f.text(0.06, 0.5, r"$\mathrm{\log \, c_{vir} (1+z)}$", fontsize=18, ha='center', va='center', rotation='vertical')
        for i in range(len(x)):
            axarr[0].errorbar(x[i],y[i],yerr=[sigy[i]],xerr=[sigx[i]],color=pl_col)
        for i in range(len(xout)):
            axarr[1].errorbar(xout[i],yout[i],yerr=[sigyout[i]],xerr=[sigxout[i]],color='k')
        plt.show()
    else:
        return xout,yout,sigxout,sigyout,uniques

def discover_repeats_all(fname_list = None, plotrepeats=False, plotdiff=False):

    if type(fname_list) is list and len(fname_list) >0:
        x,y,sigx,sigy,clusters,methods = ([],[],[],[],[],[])
        for fh in fname_list:
            x_tmp,y_tmp,sigx_tmp,sigy_tmp,clusters_tmp = startup(fname=fh)
            x = x + list(x_tmp)
            y = y + list(y_tmp)
            sigx = sigx + list(sigx_tmp)
            sigy = sigy + list(sigy_tmp)
            clusters = clusters + clusters_tmp.tolist()
            if 'X-ray' in fh:
                methods = methods + ['X-ray' for i in range(len(x_tmp))]
            elif 'WL' in fh:
                methods = methods + ['WL' for i in range(len(x_tmp))]
            elif 'SL' in fh:
                methods = methods + ['SL' for i in range(len(x_tmp))]
            elif 'WL+SL' in fh:
                methods = methods + ['WL+SL' for i in range(len(x_tmp))]
            elif 'CM' in fh:
                methods = methods + ['CM' for i in range(len(x_tmp))]
            elif 'LOSVD' in fh:
                methods = methods + ['LOSVD' for i in range(len(x_tmp))]
        assert len(methods) == len(x)
        uniques = [i for i in set(clusters)]
        counts = [clusters.count(i) for i in uniques]

        methodnames = ['X-ray','WL','SL','WL+SL','CM','LOSVD']

        xout = []
        yout = []
        sigxout = []
        sigyout = []
        methodsout = []
        for i in range(len(uniques)):
            if counts[i] == 1:
                tmp_index = clusters.index(uniques[i])
                xout.append(x[tmp_index])
                yout.append(y[tmp_index])
                sigxout.append(sigx[tmp_index])
                sigyout.append(sigy[tmp_index])
                methodsout.append(methods[tmp_index])
            else:
                tmp_indices = [j for j in range(len(clusters)) if clusters[j] == uniques[i]]
                xrep = [x[i] for i in tmp_indices]
                sigxrep = [sigx[i] for i in tmp_indices]
                xweights = [(1.0/sigxrep[i]**2) for i in range(len(xrep))]
                xnew = sum([xrep[i]*xweights[i] for i in range(len(xrep))])/sum(xweights)
                sigxnew = 1.0/np.sqrt(sum(xweights))
                xout.append(xnew)
                sigxout.append(sigxnew)
                yrep = [y[i] for i in tmp_indices]
                sigyrep = [sigy[i] for i in tmp_indices]
                yweights = [(1.0/sigyrep[i]**2) for i in range(len(yrep))]
                ynew = sum([yrep[i]*yweights[i] for i in range(len(yrep))])/sum(yweights)
                sigynew = 1.0/np.sqrt(sum(yweights))
                yout.append(ynew)
                sigyout.append(sigynew)
                tmp_methods = [methods[i] for i in tmp_indices]
                tmp_overlap = list((collections.Counter(methodnames) & collections.Counter(tmp_methods)).elements())
                if len(tmp_overlap) > 1:
                    methodsout.append('COMB')
                elif len(tmp_overlap) == 1:
                    methodsout.append(tmp_overlap[0])
                else:
                    print("No methods in tmp_overlap...")
                    return
                if plotrepeats is True:
                    plt.title("{}".format(uniques[i]))
                    plt.xlabel(r"$\mathrm{\log M_{vir}/M_{\odot}}$",fontsize=18)
                    plt.ylabel(r"$\mathrm{\log c_{vir} (1+z)}$",fontsize=18)
                    plt.errorbar(xnew,ynew,yerr=[sigynew],xerr=[sigxnew],color='cyan')
                    for i in range(len(tmp_indices)):
                        plt.errorbar(xrep[i],yrep[i],yerr=[sigyrep[i]],xerr=[sigxrep[i]],color=pl_col)
                    plt.show()
                    ipdb.set_trace()
        if plotdiff is True:
            f, axarr = plt.subplots(2, sharex=True)
            f.text(0.5, 0.01, r"$\mathrm{\log \, M_{vir}/M_{\odot}}$", fontsize=18, ha='center', va='center')
            f.text(0.06, 0.5, r"$\mathrm{\log \, c_{vir} (1+z)}$", fontsize=18, ha='center', va='center', rotation='vertical')
            for i in range(len(x)):
                axarr[0].errorbar(x[i],y[i],yerr=[sigy[i]],xerr=[sigx[i]],color=pl_col)
            for i in range(len(xout)):
                axarr[1].errorbar(xout[i],yout[i],yerr=[sigyout[i]],xerr=[sigxout[i]],color='k')
            plt.show()
        else:
            return xout,yout,sigxout,sigyout,uniques,methodsout

# Housing for a bunch of code which runs the routine for finding
# and coadding repeat measurements of clusters within each individual
# method.
def do_individual_repeat_analyses():
    # Discovering and accounting for ("co-adding") measurements
    # from the same cluster within each method, so they are not
    # over-represented in the fitting. Fitting functions in this
    # script are now accounting for repeat measurements.
    '''
    discover_repeats(method='x-ray', plotrepeats=True, plotdiff=True)
    discover_repeats(method='wl', plotrepeats=True, plotdiff=True)
    discover_repeats(method='sl', plotrepeats=True, plotdiff=True)
    discover_repeats(method='wl+sl', plotrepeats=True, plotdiff=True)
    discover_repeats(method='cm', plotrepeats=True, plotdiff=True)
    discover_repeats(method='losvd', plotrepeats=True, plotdiff=True)
    '''

# Housing for a bunch of code to run individual fitting routines on each method,
# and also makes plots of the fits (extrapolated/non-extrapolated). Also has capability to
# plot fits and uncertainties for simulations and bootstraps. Kind of sloppy at the moment,
# and it requires reading comments about what to plot (and hence uncomment).
def plot_fit_summary(extrap = False, regularsimdata = False, projectedsimdata = False, justlensing=False):
    # Plotting the fit (and uncertainty regions) over the data
    #'''
    fname = 'ObsCM_'
    if justlensing is False:
        m_xray,sigm_xray,b_xray,sigb_xray,sig_xray,(a1_xray,b1_xray,theta_xray),(xray_max,xray_min)=fit(method='xray')#,plot=True, savefig=True)
        m_sl,sigm_sl,b_sl,sigb_sl,sig_sl,(a1_sl,b1_sl,theta_sl),(sl_max,sl_min)=fit(method='sl')#,plot=True, savefig=True)
        m_cm,sigm_cm,b_cm,sigb_cm,sig_cm,(a1_cm,b1_cm,theta_cm),(cm_max,cm_min)=fit(method='cm')#,plot=True, savefig=True)
        m_losvd,sigm_losvd,b_losvd,sigb_losvd,sig_losvd,(a1_losvd,b1_losvd,theta_losvd),(losvd_max,losvd_min)=fit(method='losvd')#,plot=True, savefig=True)
        fname = fname + "All"
    else:
        fname = fname + "WLandWLSL"

    m_wl,sigm_wl,b_wl,sigb_wl,sig_wl,(a1_wl,b1_wl,theta_wl),(wl_max,wl_min)=fit(method='wl')#,plot=True, savefig=True)
    m_wlsl,sigm_wlsl,b_wlsl,sigb_wlsl,sig_wlsl,(a1_wlsl,b1_wlsl,theta_wlsl),(wlsl_max,wlsl_min)=fit(method='wl+sl')#,plot=True, savefig=True)

    ## Plotting all linear fits on the same plot
    plt.figure(figsize=(8,8),dpi=120)
    # setting range for x based upon whether or not I'm using extrapolation
    if justlensing is False:
        if extrap is False:
            xlist_xray = np.linspace(xray_min,xray_max,100)
        elif extrap is True:
            xlist_xray = np.linspace(13,17,100)
        if extrap is False:
            xlist_sl = np.linspace(sl_min,sl_max,100)
        elif extrap is True:
            xlist_sl = np.linspace(13,17,100)
        if extrap is False:
            xlist_cm = np.linspace(cm_min,cm_max,100)
        elif extrap is True:
            xlist_cm = np.linspace(13,17,100)
        if extrap is False:
            xlist_losvd = np.linspace(losvd_min,losvd_max,100)
        elif extrap is True:
            xlist_losvd = np.linspace(13,17,100)
    if extrap is False:
        fname = fname + "_NoExtrap"
        xlist_wl = np.linspace(wl_min,wl_max,100)
    elif extrap is True:
        fname = fname + "_Extrap"
        xlist_wl = np.linspace(13,17,100)
    if extrap is False:
        xlist_wlsl = np.linspace(wlsl_min,wlsl_max,100)
    elif extrap is True:
        xlist_wlsl = np.linspace(13,17,100)
    #'''
    # plotting trends and error regions for method of least-squares (w/ errors); do not plot at the same time as bootstrap section below
    #'''
    plt.plot(xlist_wl,[m_wl*i+b_wl for i in xlist_wl],color='purple',label='WL')
    plt.fill_between(xlist_wl,[(m_wl+sigm_wl)*i+(b_wl+sigb_wl+sig_wl) for i in xlist_wl],[(m_wl-sigm_wl)*i+(b_wl-sigb_wl-sig_wl) for i in xlist_wl],alpha=0.25,color='purple')
    plt.plot(xlist_wlsl,[m_wlsl*i+b_wlsl for i in xlist_wlsl],color='black',label='WL+SL')
    plt.fill_between(xlist_wlsl,[(m_wlsl+sigm_wlsl)*i+(b_wlsl+sigb_wlsl+sig_wlsl) for i in xlist_wlsl],[(m_wlsl-sigm_wlsl)*i+(b_wlsl-sigb_wlsl-sig_wlsl) for i in xlist_wlsl],alpha=0.25,color='black')
    
    if justlensing is False:
        
        plt.plot(xlist_sl,[m_sl*i+b_sl for i in xlist_sl],color='red',label='SL')
        plt.fill_between(xlist_sl,[(m_sl+sigm_sl)*i+(b_sl+sigb_sl+sig_sl) for i in xlist_sl],[(m_sl-sigm_sl)*i+(b_sl-sigb_sl-sig_sl) for i in xlist_sl],alpha=0.25,color='red')
        plt.plot(xlist_xray,[m_xray*i+b_xray for i in xlist_xray],color='green',label='X-ray')
        plt.fill_between(xlist_xray,[(m_xray+sigm_xray)*i+(b_xray+sigb_xray+sig_xray) for i in xlist_xray],[(m_xray-sigm_xray)*i+(b_xray-sigb_xray-sig_xray) for i in xlist_xray],alpha=0.25,color='green')
        plt.plot(xlist_cm,[m_cm*i+b_cm for i in xlist_cm],color='blue',label='CM')
        plt.fill_between(xlist_cm,[(m_cm+sigm_cm)*i+(b_cm+sigb_cm+sig_cm) for i in xlist_cm],[(m_cm-sigm_cm)*i+(b_cm-sigb_cm-sig_cm) for i in xlist_cm],alpha=0.25,color='blue')
        plt.plot(xlist_losvd,[m_losvd*i+b_losvd for i in xlist_losvd],color='orange',label='LOSVD')
        plt.fill_between(xlist_losvd,[(m_losvd+sigm_losvd)*i+(b_losvd+sigb_losvd+sig_losvd) for i in xlist_losvd],[(m_losvd-sigm_losvd)*i+(b_losvd-sigb_losvd-sig_losvd) for i in xlist_losvd],alpha=0.25,color='orange')
    #'''
    # plotting bootstrap fits and error regions (w/ errors), manually
    '''
    plt.plot(xlist_xray,[-0.1627*i+3.3331 for i in xlist_xray],color='green',label='X-ray')
    plt.fill_between(xlist_xray,[(-0.1627+0.0284)*i+(3.3348+0.4274+0.1994) for i in xlist_xray],[(-0.1627-0.0284)*i+(3.3348-0.4274-0.1994) for i in xlist_xray],alpha=0.25,color='green')
    plt.plot(xlist_wl,[-0.4917*i+8.2976 for i in xlist_wl],color='purple',label='WL')
    plt.fill_between(xlist_wl,[(-0.4917+0.0953)*i+(8.2976+1.4125+0.1716) for i in xlist_wl],[(-0.4917-0.0953)*i+(8.2976-1.4125-0.1716) for i in xlist_wl],alpha=0.25,color='purple')
    plt.plot(xlist_sl,[-0.0003*i+0.9736 for i in xlist_sl],color='red',label='SL')
    plt.fill_between(xlist_sl,[(-0.0003+0.1693)*i+(0.9736+2.5356+0.1188) for i in xlist_sl],[(-0.0003-0.1693)*i+(0.9736-2.5356-0.1188) for i in xlist_sl],alpha=0.25,color='red')
    plt.plot(xlist_wlsl,[-0.6043*i+9.9818 for i in xlist_wlsl],color='black',label='WL+SL')
    plt.fill_between(xlist_wlsl,[(-0.6043+0.1325)*i+(9.9818+1.9602+0.1411) for i in xlist_wlsl],[(-0.6043-0.1325)*i+(9.9818-1.9602-0.1411) for i in xlist_wlsl],alpha=0.25,color='black')
    plt.plot(xlist_cm,[0.2187*i+-2.2660 for i in xlist_cm],color='blue',label='CM')
    plt.fill_between(xlist_cm,[(0.2187+0.1196)*i+(-2.2660+1.7296+0.2155) for i in xlist_cm],[(0.2187-0.1196)*i+(-2.2660-1.7296-0.2155) for i in xlist_cm],alpha=0.25,color='blue')
    plt.plot(xlist_losvd,[0.0379*i+0.2463 for i in xlist_losvd],color='orange',label='LOSVD')
    plt.fill_between(xlist_losvd,[(0.0379+0.1225)*i+(0.2463+1.8035+0.1797) for i in xlist_losvd],[(0.0379-0.1225)*i+(0.2463-1.8035-0.1797) for i in xlist_losvd],alpha=0.25,color='orange')
    '''
    ## plotting sim fits over-top of data
    #m_sim_los_half,b_sim_los_half,sig_sim_los_half=fit_sims(los_angle='aligned',scale='half')
    #m_sim_los_full,b_sim_los_full,sig_sim_los_full=fit_sims(los_angle='aligned',scale='full')
    #m_sim_intrinsic,b_sim_intrinsic,sig_sim_intrinsic=fit_sims()
    #plt.plot(np.linspace(13,17,100),[m_sim_intrinsic*i+b_sim_intrinsic for i in np.linspace(13,17,100)],color='cyan',label='Intrinsic Simulations')
    #plt.fill_between(np.linspace(13,17,100),[(m_sim_intrinsic)*i+(b_sim_intrinsic+sig_sim_intrinsic) for i in np.linspace(13,17,100)],[(m_sim_intrinsic)*i+(b_sim_intrinsic-sig_sim_intrinsic) for i in np.linspace(13,17,100)],alpha=0.25,color='cyan')

    ## Plotting sim data over-top of observations
    #'''
    if regularsimdata is True:
        l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()
        plt.errorbar(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),yerr=np.log10(np.std(l_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,label='Groener and Goldberg (2014)',zorder=20,alpha=1.0)
        plt.scatter(np.log10(np.average(l_masses)),np.log10(np.average(l_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k',alpha=1.0)
        plt.errorbar(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),yerr=np.log10(np.std(m_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20,alpha=1.0)
        plt.scatter(np.log10(np.average(m_masses)),np.log10(np.average(m_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k',alpha=1.0)
        plt.errorbar(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),yerr=np.log10(np.std(h_concs)),
                     fmt='*',color='magenta',capsize=10,capthick=3,elinewidth=8,zorder=20,alpha=1.0)
        plt.scatter(np.log10(np.average(h_masses)),np.log10(np.average(h_concs)),
                    marker='*',s=500,zorder=21,color='magenta',edgecolor='k',alpha=1.0)
        fname = fname + "andSims"
    #'''
    # projecting halos
    if projectedsimdata is True:
        l_concs_p = [conc_finder_pro(l_concs[i],[0],0.5) for i in range(len(l_concs))]
        m_concs_p = [conc_finder_pro(m_concs[i],[0],0.5) for i in range(len(m_concs))]
        h_concs_p = [conc_finder_pro(h_concs[i],[0],0.5) for i in range(len(h_concs))]
        plt.errorbar(np.log10(np.average(l_masses)),np.log10(np.average(l_concs_p)),yerr=np.log10(np.std(l_concs_p)),
                     fmt='*',color='cyan',capsize=10,capthick=3,elinewidth=8,zorder=19,alpha=1.0)
        plt.scatter(np.log10(np.average(l_masses)),np.log10(np.average(l_concs_p)),
                    marker='*',s=500,zorder=23,color='cyan',edgecolor='k',alpha=1.0)
        plt.errorbar(np.log10(np.average(m_masses)),np.log10(np.average(m_concs_p)),yerr=np.log10(np.std(m_concs_p)),
                     fmt='*',color='cyan',capsize=10,capthick=3,elinewidth=8,zorder=19,alpha=1.0)
        plt.scatter(np.log10(np.average(m_masses)),np.log10(np.average(m_concs_p)),
                    marker='*',s=500,zorder=23,color='cyan',edgecolor='k',alpha=1.0)
        plt.errorbar(np.log10(np.average(h_masses)),np.log10(np.average(h_concs_p)),yerr=np.log10(np.std(h_concs_p)),
                     fmt='*',color='cyan',capsize=10,capthick=3,elinewidth=8,zorder=19,alpha=1.0)
        plt.scatter(np.log10(np.average(h_masses)),np.log10(np.average(h_concs_p)),
                    marker='*',s=500,zorder=23,color='cyan',edgecolor='k',alpha=1.0)
        fname = fname + "andPSims.png"
    else:
        fname = fname + '.png'

    plt.legend(loc=0,numpoints=1,frameon=False,fontsize=12)
    plt.xlabel(r'$\mathrm{\log \, M_{vir}/M_{\odot}}$',fontsize=25)
    plt.ylabel(r'$\mathrm{\log \, c_{vir} (1+z) }$',fontsize=25)
    
    plt.savefig(fname)
    
    
    # Plotting all of the error ellipses on the same plot
    # Not that informative of a plot...
    '''
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)

    ell_xray=Ellipse([m_xray,b_xray],width=2*a1_xray,height=2*b1_xray,angle=theta_xray*57.3)
    ell2_xray=Ellipse([m_xray,b_xray],width=a1_xray,height=b1_xray,angle=theta_xray*57.3)
    ax.add_artist(ell_xray)
    ell_xray.set_facecolor('white')
    ell_xray.set_edgecolor('green')
    ax.add_artist(ell2_xray)
    ell2_xray.set_facecolor('white')
    ell2_xray.set_edgecolor('green')
    ax.plot(m_xray,b_xray,'+',color='green')

    ell_wl=Ellipse([m_wl,b_wl],width=2*a1_wl,height=2*b1_wl,angle=theta_wl*57.3)
    ell2_wl=Ellipse([m_wl,b_wl],width=a1_wl,height=b1_wl,angle=theta_wl*57.3)
    ax.add_artist(ell_wl)
    ell_wl.set_facecolor('white')
    ell_wl.set_edgecolor('purple')
    ax.add_artist(ell2_wl)
    ell2_wl.set_facecolor('white')
    ell2_wl.set_edgecolor('purple')
    ax.plot(m_wl,b_wl,'+',color='purple')

    ell_sl=Ellipse([m_sl,b_sl],width=2*a1_sl,height=2*b1_sl,angle=theta_sl*57.3)
    ell2_sl=Ellipse([m_sl,b_sl],width=a1_sl,height=b1_sl,angle=theta_sl*57.3)
    ax.add_artist(ell_sl)
    ell_sl.set_facecolor('white')
    ell_sl.set_edgecolor('red')
    ax.add_artist(ell2_sl)
    ell2_sl.set_facecolor('white')
    ell2_sl.set_edgecolor('red')
    ax.plot(m_sl,b_sl,'+',color='red')

    ell_wlsl=Ellipse([m_wlsl,b_wlsl],width=2*a1_wlsl,height=2*b1_wlsl,angle=theta_wlsl*57.3)
    ell2_wlsl=Ellipse([m_wlsl,b_wlsl],width=a1_wlsl,height=b1_wlsl,angle=theta_wlsl*57.3)
    ax.add_artist(ell_wlsl)
    ell_wlsl.set_facecolor('white')
    ell_wlsl.set_edgecolor('black')
    ax.add_artist(ell2_wlsl)
    ell2_wlsl.set_facecolor('white')
    ell2_wlsl.set_edgecolor('black')
    ax.plot(m_wlsl,b_wlsl,'+',color='black')

    ell_cm=Ellipse([m_cm,b_cm],width=2*a1_cm,height=2*b1_cm,angle=theta_cm*57.3)
    ell2_cm=Ellipse([m_cm,b_cm],width=a1_cm,height=b1_cm,angle=theta_cm*57.3)
    ax.add_artist(ell_cm)
    ell_cm.set_facecolor('white')
    ell_cm.set_edgecolor('blue')
    ax.add_artist(ell2_cm)
    ell2_cm.set_facecolor('white')
    ell2_cm.set_edgecolor('blue')
    ax.plot(m_cm,b_cm,'+',color='blue')

    ell_losvd=Ellipse([m_losvd,b_losvd],width=2*a1_losvd,height=2*b1_losvd,angle=theta_losvd*57.3)
    ell2_losvd=Ellipse([m_losvd,b_losvd],width=a1_losvd,height=b1_losvd,angle=theta_losvd*57.3)
    ax.add_artist(ell_losvd)
    ell_losvd.set_facecolor('white')
    ell_losvd.set_edgecolor('orange')
    ax.add_artist(ell2_losvd)
    ell2_losvd.set_facecolor('white')
    ell2_losvd.set_edgecolor('orange')
    ax.plot(m_losvd,b_losvd,'+',color='orange')

    plt.xlabel('m')
    plt.ylabel('b')

    plt.show()
    '''

# Housing for a bunch of code to run bootstrap analyses for each method (with and without errors)
def boostrap_summary():
    # Performing bootstrap
    #tmp_path = '/Users/groenera/Desktop/Dropbox/Private/Research/GroupMeetings/Meeting#60/' # on OS X
    tmp_path = '/home/groenera/Desktop/Dropbox/Private/Research/GroupMeetings/Meeting#60/' # on Ubuntu
    #without error first
    '''
    do_bootstrap(method='X-ray',witherrors=False, savepath=tmp_path)
    do_bootstrap(method='CM',witherrors=False, savepath=tmp_path)
    do_bootstrap(method='WL',witherrors=False, savepath=tmp_path)
    do_bootstrap(method='SL',witherrors=False, savepath=tmp_path)
    do_bootstrap(method='WL+SL',witherrors=False, savepath=tmp_path)
    do_bootstrap(method='LOSVD',witherrors=False, savepath=tmp_path)
    '''
    #with error next
    #'''
    #do_bootstrap(method='X-ray',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='CM',witherrors=True,savepath=tmp_path,nsamples=200)
    do_bootstrap(method='WL',witherrors=True,savepath=tmp_path,nsamples=100)
    do_bootstrap(method='SL',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='WL+SL',witherrors=True,savepath=tmp_path,nsamples=200)
    #do_bootstrap(method='LOSVD',witherrors=True,savepath=tmp_path,nsamples=100)
    #'''

def plot_sample_summary(plotrepeats=True, savefigure=True, witherrors=True):
    if plotrepeats is False:
        x_xray,y_xray,sigx_xray,sigy_xray,cl_xray = discover_repeats(method='X-ray')
        x_wl,y_wl,sigx_wl,sigy_wl,cl_wl = discover_repeats(method='WL')
        x_sl,y_sl,sigx_sl,sigy_sl,cl_sl = discover_repeats(method='SL')
        x_wlsl,y_wlsl,sigx_wlsl,sigy_wlsl,cl_wlsl = discover_repeats(method='WL+SL')
        x_cm,y_cm,sigx_cm,sigy_cm,cl_cm = discover_repeats(method='CM')
        x_losvd,y_losvd,sigx_losvd,sigy_losvd,cl_losvd = discover_repeats(method='LOSVD')
    elif plotrepeats is True:
        x_xray,y_xray,sigx_xray,sigy_xray,cl_xray = startup(fname='X-ray_data.txt')
        x_wl,y_wl,sigx_wl,sigy_wl,cl_wl = startup(fname='WL_data.txt')
        x_sl,y_sl,sigx_sl,sigy_sl,cl_sl = startup(fname='SL_data.txt')
        x_wlsl,y_wlsl,sigx_wlsl,sigy_wlsl,cl_wlsl = startup(fname='WL+SL_data.txt')
        x_cm,y_cm,sigx_cm,sigy_cm,cl_cm = startup(fname='CM_data.txt')
        x_losvd,y_losvd,sigx_losvd,sigy_losvd,cl_losvd = startup(fname='LOSVD_data.txt')
    title = "PLOT SAMPLE SUMMARY"
    print(len(title)*'*')
    print(title)
    print(len(title)*'*')
    print("Number of unique clusters for X-ray: {}".format(len(x_xray)))
    print("Number of unique clusters for WL: {}".format(len(x_wl)))
    print("Number of unique clusters for SL: {}".format(len(x_sl)))
    print("Number of unique clusters for WL+SL: {}".format(len(x_wlsl)))
    print("Number of unique clusters for CM: {}".format(len(x_cm)))
    print("Number of unique clusters for LOSVD: {}".format(len(x_losvd)))
    print("Number of total unique clusters: {}".format(len(x_xray)+len(x_wl)+len(x_sl)+len(x_wlsl)+len(x_cm)+len(x_losvd)))
    plt.figure(1,figsize=(7.5,7.5))
    for i in range(len(x_xray)):
        if witherrors is True:
            plt.errorbar(x_xray[i],y_xray[i],yerr=sigy_xray[i],xerr=sigx_xray[i],color='green')
        elif witherrors is False:
            plt.scatter(x_xray[i],y_xray[i],marker='*',color='green')
    for i in range(len(x_wl)):
        if witherrors is True:
            plt.errorbar(x_wl[i],y_wl[i],yerr=sigy_wl[i],xerr=sigx_wl[i],color='purple',zorder=2)
        elif witherrors is False:
            plt.scatter(x_wl[i],y_wl[i],color='purple',marker='d',zorder=2)
    for i in range(len(x_sl)):
        if witherrors is True:
            plt.errorbar(x_sl[i],y_sl[i],yerr=sigy_sl[i],xerr=sigx_sl[i],color='red')
        elif witherrors is False:
            plt.scatter(x_sl[i],y_sl[i],marker='s',color='red')
    for i in range(len(x_wlsl)):
        if witherrors is True:
            plt.errorbar(x_wlsl[i],y_wlsl[i],yerr=sigy_wlsl[i],xerr=sigx_wlsl[i],color='black')
        elif witherrors is False:
            plt.scatter(x_wlsl[i],y_wlsl[i],marker='o',color='black')
    for i in range(len(x_cm)):
        if witherrors is True:
            plt.errorbar(x_cm[i],y_cm[i],yerr=sigy_cm[i],xerr=sigx_cm[i],color='blue')
        elif witherrors is False:
            plt.scatter(x_cm[i],y_cm[i],marker='x',color='blue')
    for i in range(len(x_losvd)):
        if witherrors is True:
            plt.errorbar(x_losvd[i],y_losvd[i],yerr=sigy_losvd[i],xerr=sigx_losvd[i],color='orange')
        elif witherrors is False:
            plt.scatter(x_losvd[i],y_losvd[i],marker='^',color='orange')
    if witherrors is True:
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='green',label='X-ray')
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='purple',label='WL')
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='red',label='SL')
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='black',label='WL+SL')
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='blue',label='CM')
        plt.errorbar(1e16,1e16,yerr=1e16,xerr=1e16,color='orange',label='LOSVD')
    elif witherrors is False:
        plt.scatter(1e16,1e16,color='green',marker='*',label='X-ray')
        plt.scatter(1e16,1e16,color='purple',marker='d',label='WL')
        plt.scatter(1e16,1e16,color='red',marker='s',label='SL')
        plt.scatter(1e16,1e16,color='black',marker='o',label='WL+SL')
        plt.scatter(1e16,1e16,color='blue',marker='x',label='CM')
        plt.scatter(1e16,1e16,color='orange',marker='^',label='LOSVD')
    plt.xlim(13.0,17.0)
    plt.ylim(-0.5,2.0)
    plt.legend(loc=0,numpoints=1,scatterpoints=1,fontsize=12,frameon=True)
    plt.xlabel(r'$\mathrm{\log\,{ M_{vir}/M_{\odot}}}$',fontsize=18)
    plt.ylabel(r'$\mathrm{\log\,{ \, c_{vir} \, (1+z) }}$',fontsize=18)
    if savefigure is True:
        plt.savefig("CMRelation_FullSample_Symmetrized.png")
    plt.tight_layout()
    plt.show()

    
if __name__ == "__main__":

    # Doing fits to individual fits to each method sub-population
    # If plot=True and savefig=True, figures are saved in the current directory
    #fit(method='X-ray', plot=True, savefig=True)
    #fit(method='WL', plot=True, savefig=True)
    #for strong lensing, use bootstrap plot instead
    #fit(method='SL', plot=True, savefig=True)
    '''
    m_list, b_list, sig = fit_bootstrap(method='sl', witherrors=True, nsamples=100)
    m_ave,m_std = (np.average(m_list),np.std(m_list))
    b_ave,b_std = (np.average(b_list),np.std(b_list))
    x,y,sigx,sigy,cl=discover_repeats(method='sl')
    plt.figure(figsize=(8,8))
    plt.title('SL')
    plt.errorbar(x,y,xerr=sigx,yerr=sigy,fmt='o',color='red')
    x0 = np.array([13.0,17.5])
    y0=m_ave*x0+b_ave
    plt.plot(x0,y0,linewidth=3,color='red')
    plt.fill_between(x0,[(m_ave+m_std)*i+(b_ave+b_std+sig) for i in x0],[(m_ave-m_std)*i+(b_ave-b_std-sig) for i in x0],alpha=0.25,color='red')
    plt.xlabel(r'$\mathrm{\log{ M_{vir}/M_{\odot}}}$',fontsize=18)
    plt.ylabel(r'$\mathrm{\log{ \, c_{vir} \, (1+z) }}$',fontsize=18)
    plt.xlim(13.0,17.5)
    plt.ylim(-1.0,2.5)
    plt.show()
    '''
    #fit(method='WL+SL', plot=True, savefig=True)
    #fit(method='CM', plot=True, savefig=True)
    #fit(method='LOSVD', plot=True, savefig=True)


    # Making plot of fit to ALL data (with data plotted, too),
    # with sims overlaid on top.
    #fit_all(plot=True, savefig=True, plotwitherrorbars=True)

    # Making plot of fits to WL and WL+SL individually
    # (no individual data points plotted here), with sims overlaid
    # on top
    #plot_fit_summary(extrap=False, regularsimdata=True, projectedsimdata=True, justlensing=True)
    #plot_fit_summary(extrap=False, regularsimdata=False, projectedsimdata=False, justlensing=False)
    
    # Making plot of the full sample (masses/concs)
    #plot_sample_summary(witherrors=False)

    # Doing full bootstrap analysis on all methods (for each method individually)
    boostrap_summary()

    # Doing full bootstrap analysis on the full sample at once
    #m2_list, b2_list, sig = fit_bootstrap_allmethods(witherrors=True, nsamples=100)
    
