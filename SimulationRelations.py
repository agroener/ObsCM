import matplotlib.pyplot as plt
import numpy as np

# Simulation relation functions
def bullock(Mvir,z,Mstar=1.3e13):
    '''
    Title: Profiles of dark haloes: evolution, scatter, and environment
    ADS Link: http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:astro-ph/9908159
    Comments: scatter of 0.18 in log(c_vir)
    '''
    return (9.0/(1+z)) * (Mvir/Mstar)**(-0.13)
def hannawi(Mvir,z,Mstar=1.3e13):
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return (12.3/(1+z)) * (Mvir/Mstar)**(-0.13)    
def prada():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def zhao():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def ludlow():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def dutton():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def battacharya():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def correa():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def diemer():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return
def vandenbosch():
    '''
    Title: 
    ADS Link: 
    Comments: 
    '''
    return

# Function for finding Groener and Goldberg data
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

    return l_concs,l_masses,m_concs,m_masses,h_concs,h_masses


if __name__ == "__main__":
    
    # Get simulation results from Groener & Goldberg 2014, first
    l_concs,l_masses,m_concs,m_masses,h_concs,h_masses = startup_sims()

    # Call various simulation relations here
