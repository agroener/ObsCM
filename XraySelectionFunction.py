import numpy as np
import matplotlib.pyplot as plt
import random

#################
### Functions ###
#################

# Only Flat Cosmologies
def Ez(z,Omega_m_0 = 0.3, Omega_L_0 = 0.7):
    return np.sqrt(Omega_m_0*(1+z)**3 + Omega_L_0)

