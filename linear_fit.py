from numpy import *
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.misc import derivative
import numpy as np
import ipdb

# Functions

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
        x,y,sigx,sigy=loadtxt(fname,delimiter=",",usecols=(0,1,2,3),unpack=True,comments='#')
        return x,y,sigx,sigy
    else:
        print("Filename must be specified...")
        return

def startup_bootstrap(fname=None):
    # first load in all of the data
    if fname is not None:
        x,y,sigx,sigy=loadtxt(fname,delimiter=",",usecols=(0,1,2,3),unpack=True,comments='#')
    else:
        print("Filename must be specified...")
        return
    # now sample with replacement
    x_bs = np.random.choice(x,replace=True,size=len(x))
    bs_indices = [np.where(x == x_bs[i])[0][0] for i in range(len(x_bs))]
    y_bs = np.array([y[i] for i in bs_indices])
    sigx_bs = np.array([sigx[i] for i in bs_indices])
    sigy_bs = np.array([sigy[i] for i in bs_indices])
    return x_bs,y_bs,sigx_bs,sigy_bs
        
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
    
    x,y,sigx,sigy = startup(fh=filename)

    print("Number of Measurements Used: {}".format(len(x)))
    
    # Fitting as if there are no uncertainties
    N=len(x)
    Sxy=sum(x*y)
    Sx=sum(x)
    Sy=sum(y)
    Sxx=sum(x*x)
    m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
    b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
    sig=sqrt(std(y-(m1*x+b1))**2-mean(sigy**2))
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
    sigm=1./sqrt(F11)
    sigb=1./sqrt(F22)

    print("Linear model (with uncertainties): m: {} +/- {}, b: {} +/- {}".format(m2,sigm,b2,sigb))

    
    # Plotting the best fit line over the data
    if plot is True:
        plt.figure(figsize=(8,8))
        plt.title(method_title)
        plt.errorbar(x,y,xerr=sigx,yerr=sigy,fmt='o',color=pl_col)
        x0 = array([13.0,17.5])
        y0=m2*x0+b2
        plt.plot(x0,y0,linewidth=3,color=pl_col)
        plt.fill_between(x0,[(m2+sigm)*i+(b2+sigb) for i in x0],[(m2-sigm)*i+(b2-sigb) for i in x0],alpha=0.25,color=pl_col)
        plt.xlabel(r'$\mathrm{\log{ M_{vir}/M_{\odot}}}$',fontsize=18)
        plt.ylabel(r'$\mathrm{\log{ \, c_{vir} \cdot (1+z) }}$',fontsize=18)
        plt.xlim(13.0,17.5)
        plt.ylim(-1.0,2.5)
        if savefig is True:
            plt.savefig('{}_linearmodel_witherror.png'.format(filename.split('_')[0]))

    # Calculating and plotting error ellipse
    detF=F11*F22-F12*F12
    Qxx=F22/detF
    Qyy=F11/detF
    Qxy=-F12/detF
    
    theta=arctan(2*Qxy/(Qxx-Qyy+1e-9))/2.
    a1=sqrt(Qxx*pow(cos(theta),2)+Qyy*pow(sin(theta),2)+2*Qxy*sin(theta)*cos(theta))
    b1=sqrt(Qxx*pow(sin(theta),2)+Qyy*pow(cos(theta),2)-2*Qxy*sin(theta)*cos(theta))    

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

    return m2,sigm,b2,sigb,sig,(a1,b1,theta)

def fit_bootstrap(method='X-ray', witherrors=True, nsamples=100):
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
            x_bs,y_bs,sigx_bs,sigy_bs = startup_bootstrap(fname=filename)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=sqrt(std(y_bs-(m1*x_bs+b1))**2-mean(sigy_bs**2))
            m1_list.append(m1)
            b1_list.append(b1)
            sig_list.append(sig)
        return m1_list, b1_list, sig_list

    if witherrors is True:
        m2_list, b2_list = ([],[])
        for i in range(nsamples):
            x_bs,y_bs,sigx_bs,sigy_bs = startup_bootstrap(fname=filename)
            N=len(x_bs)
            Sxy=sum(x_bs*y_bs)
            Sx=sum(x_bs)
            Sy=sum(y_bs)
            Sxx=sum(x_bs*x_bs)
            m1=(N*Sxy-Sx*Sy)/(N*Sxx-Sx*Sx)
            b1=(-Sx*Sxy+Sxx*Sy)/(N*Sxx-Sx*Sx)
            sig=sqrt(abs(std(y_bs-(m1*x_bs+b1))**2-mean(sigy_bs**2))) #abs needed for cases where second term is larger
            m2,b2 = steepest_decent(x_bs,y_bs,sigx_bs,sigy_bs,sig,m1,b1,N,alpha=0.75,tol=1.e-6)
            m2_list.append(m2)
            b2_list.append(b2)
        return m2_list, b2_list

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
        m, b = fit_bootstrap(method=method, witherrors=witherrors, nsamples=nsamples)
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

    
if __name__ == "__main__":

    # Plotting the fit (and uncertainty regions) over the data
    '''
    m_xray,sigm_xray,b_xray,sigb_xray,sig_xray,(a1_xray,b1_xray,theta_xray)=fit(method='xray')#,plot=True, savefig=True)
    m_wl,sigm_wl,b_wl,sigb_wl,sig_wl,(a1_wl,b1_wl,theta_wl)=fit(method='wl')#,plot=True, savefig=True)
    m_sl,sigm_sl,b_sl,sigb_sl,sig_sl,(a1_sl,b1_sl,theta_sl)=fit(method='sl')#,plot=True, savefig=True)
    m_wlsl,sigm_wlsl,b_wlsl,sigb_wlsl,sig_wlsl,(a1_wlsl,b1_wlsl,theta_wlsl)=fit(method='wl+sl')#,plot=True, savefig=True)
    m_cm,sigm_cm,b_cm,sigb_cm,sig_cm,(a1_cm,b1_cm,theta_cm)=fit(method='cm')#,plot=True, savefig=True)
    m_losvd,sigm_losvd,b_losvd,sigb_losvd,sig_losvd,(a1_losvd,b1_losvd,theta_losvd)=fit(method='losvd')#,plot=True, savefig=True)

    # Plotting all linear fits on the same plot
    xlist = np.linspace(13,17,300)
    plt.figure(figsize=(8,8))
    plt.plot(xlist,[m_xray*i+b_xray for i in xlist],color='green',label='X-ray')
    plt.fill_between(xlist,[(m_xray+sigm_xray)*i+(b_xray+sigb_xray+sig_xray) for i in xlist],[(m_xray-sigm_xray)*i+(b_xray-sigb_xray-sig_xray) for i in xlist],alpha=0.25,color='green')
    plt.plot(xlist,[m_wl*i+b_wl for i in xlist],color='purple',label='WL')
    plt.fill_between(xlist,[(m_wl+sigm_wl)*i+(b_wl+sigb_wl+sig_wl) for i in xlist],[(m_wl-sigm_wl)*i+(b_wl-sigb_wl-sig_wl) for i in xlist],alpha=0.25,color='purple')
    plt.plot(xlist,[m_sl*i+b_sl for i in xlist],color='red',label='SL')
    plt.fill_between(xlist,[(m_sl+sigm_sl)*i+(b_sl+sigb_sl+sig_sl) for i in xlist],[(m_sl-sigm_sl)*i+(b_sl-sigb_sl-sig_sl) for i in xlist],alpha=0.25,color='red')
    plt.plot(xlist,[m_wlsl*i+b_wlsl for i in xlist],color='black',label='WL+SL')
    plt.fill_between(xlist,[(m_wlsl+sigm_wlsl)*i+(b_wlsl+sigb_wlsl+sig_wlsl) for i in xlist],[(m_wlsl-sigm_wlsl)*i+(b_wlsl-sigb_wlsl-sig_wlsl) for i in xlist],alpha=0.25,color='black')
    plt.plot(xlist,[m_cm*i+b_cm for i in xlist],color='blue',label='CM')
    plt.fill_between(xlist,[(m_cm+sigm_cm)*i+(b_cm+sigb_cm+sig_cm) for i in xlist],[(m_cm-sigm_cm)*i+(b_cm-sigb_cm-sig_cm) for i in xlist],alpha=0.25,color='blue')
    plt.plot(xlist,[m_losvd*i+b_losvd for i in xlist],color='orange',label='LOSVD')
    plt.fill_between(xlist,[(m_losvd+sigm_losvd)*i+(b_losvd+sigb_losvd+sig_losvd) for i in xlist],[(m_losvd-sigm_losvd)*i+(b_losvd-sigb_losvd-sig_losvd) for i in xlist],alpha=0.25,color='orange')
    plt.legend(loc=0)
    plt.xlabel(r'$\mathrm{\log \, M_{vir}/M_{\odot}}$',fontsize=20)
    plt.ylabel(r'$\mathrm{\log \, c_{vir} (1+z) }$',fontsize=20)
    plt.savefig('ObsCM_LinearModel.png')
    '''
    
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

    # Performing bootstrap
    tmp_path = '/Users/groenera/Desktop/Dropbox/Private/Research/GroupMeetings/Meeting#59/'
    #without error first
    #do_bootstrap(method='X-ray',witherrors=False, savepath=tmp_path)
    #do_bootstrap(method='CM',witherrors=False, savepath=tmp_path)
    #do_bootstrap(method='WL',witherrors=False, savepath=tmp_path)
    #do_bootstrap(method='SL',witherrors=False, savepath=tmp_path)
    #do_bootstrap(method='WL+SL',witherrors=False, savepath=tmp_path)
    #do_bootstrap(method='LOSVD',witherrors=False, savepath=tmp_path)
    #with error next
    #do_bootstrap(method='X-ray',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='CM',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='WL',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='SL',witherrors=True,savepath=tmp_path,nsamples=100)
    do_bootstrap(method='WL+SL',witherrors=True,savepath=tmp_path,nsamples=100)
    #do_bootstrap(method='LOSVD',witherrors=True,savepath=tmp_path,nsamples=100)
