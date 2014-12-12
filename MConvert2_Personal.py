from __future__ import division # float division by default
import numpy as np
import ipdb
import MConvert_Personal as mc

def MConvert_SE14(Mold,DeltaOld,DeltaNew,ConcentOld,redshift,method,Omega_m_0_new,Omega_L_0_new,Omega_m_0_old=0.3,Omega_L_0_old=0.7):
    Mnew = mc.Mconvert(Mold,DeltaOld,DeltaNew,ConcentOld)
    Cnew = mc.Cconvert(Mold,DeltaOld,DeltaNew,ConcentOld)
    if method is None:
        return Mnew,Cnew
    elif method in ['WL','SL','WL+SL']:
        ipdb.set_trace()
        


if __name__ == "__main__":
    MConvert_SE14(10.5,200,500,7.5,0.183,'WL',0.32,0.68)
