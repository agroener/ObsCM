import numpy as np
import matplotlib.pyplot as plt

def writedata(mvir_norm,mvir_p_norm,cvir_norm,cvir_p_norm,methods_norm,z_norm,cl_norm,method=None,plot=False,ref=None,refs_norm=None):

    if method in ['WL','wl']:
        m = 'WL'
    elif method in ['X-ray','x-ray','xray']:
        m = 'X-ray'
    elif method in ['SL','sl']:
        m = 'SL'
    elif method in ['WL+SL','wl+sl']:
        m = 'WL+SL'
    elif method in ['CM','cm']:
        m = 'CM'
    elif method in ['LOSVD','losvd']:
        m = 'LOSVD'
    else:
        print("Method is ill-defined...")
        return
    
    if ref is None:
        x = [np.log10(mvir_norm[i]*1e14) for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0]
        y = [np.log10(cvir_norm[i]*(1+z_norm[i])) for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0]
        z = [z_norm[i] for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0]
        sigx = [0.434*(mvir_p_norm[i]/mvir_norm[i]) for i in range(len(mvir_norm))
                if methods_norm[i] == m
                and mvir_norm[i]-mvir_p_norm[i]>0
                and cvir_norm[i]-cvir_p_norm[i]>0]
        sigy = [0.434*(cvir_p_norm[i]/cvir_norm[i]) for i in range(len(mvir_norm))
                if methods_norm[i] == m
                and mvir_norm[i]-mvir_p_norm[i]>0
                and cvir_norm[i]-cvir_p_norm[i]>0]
        cl = [cl_norm[i] for i in range(len(mvir_norm))
              if methods_norm[i] == m
              and mvir_norm[i]-mvir_p_norm[i]>0
              and cvir_norm[i]-cvir_p_norm[i]>0]
    else:
        x = [np.log10(mvir_norm[i]*1e14) for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0
             and refs_norm[i] == ref]
        y = [np.log10(cvir_norm[i]*(1+z_norm[i])) for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0
             and refs_norm[i] == ref]
        z = [z_norm[i] for i in range(len(mvir_norm))
             if methods_norm[i] == m
             and mvir_norm[i]-mvir_p_norm[i]>0
             and cvir_norm[i]-cvir_p_norm[i]>0
             and refs_norm[i] == ref]
        sigx = [0.434*(mvir_p_norm[i]/mvir_norm[i]) for i in range(len(mvir_norm))
                if methods_norm[i] == m
                and mvir_norm[i]-mvir_p_norm[i]>0
                and cvir_norm[i]-cvir_p_norm[i]>0
                and refs_norm[i] == ref]
        sigy = [0.434*(cvir_p_norm[i]/cvir_norm[i]) for i in range(len(mvir_norm))
                if methods_norm[i] == m
                and mvir_norm[i]-mvir_p_norm[i]>0
                and cvir_norm[i]-cvir_p_norm[i]>0
                and refs_norm[i] == ref]
        cl = [cl_norm[i] for i in range(len(mvir_norm))
              if methods_norm[i] == m
              and mvir_norm[i]-mvir_p_norm[i]>0
              and cvir_norm[i]-cvir_p_norm[i]>0
              and refs_norm[i] == ref]
        m = ref + "_" + m
        
    if plot is True:
        for i in range(len(x)):
            plt.errorbar(x[i],y[i],yerr=[sigx[i]],xerr=[sigy[i]],color='blue')
        plt.show()

    FH = open('{}_data.txt'.format(m), 'w')
    for i in range(len(x)):
        line = "{},{},{},{},{}".format(cl[i],x[i],y[i],sigx[i],sigy[i])
        FH.write(line+"\n")
    FH.close()
