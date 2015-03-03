from xlrd import open_workbook

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

import sys
import copy

from astropy.table import Table
from astropy.io import ascii

import ipdb

def startup():
    # Opening Excel File
    wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')
    # First sheet
    sheet_names = wb.sheet_names()
    sheet1 = wb.sheet_by_name(sheet_names[0])
    # Get headers and data
    headers = []
    clusters = sheet1.col_values(0)
    headers.append(clusters.pop(0))
    redshift = sheet1.col_values(1)
    headers.append(redshift.pop(0))
    methods = sheet1.col_values(2)
    headers.append(methods.pop(0))
    c200 = sheet1.col_values(3)
    headers.append(c200.pop(0))
    c200_plus = sheet1.col_values(4)
    headers.append(c200_plus.pop(0))
    c200_minus = sheet1.col_values(5)
    headers.append(c200_minus.pop(0))
    m200 = sheet1.col_values(6)
    headers.append(m200.pop(0))
    m200_plus = sheet1.col_values(7)
    headers.append(m200_plus.pop(0))
    m200_minus = sheet1.col_values(8)
    headers.append(m200_minus.pop(0))
    cvir = sheet1.col_values(9)
    headers.append(cvir.pop(0))
    cvir_plus = sheet1.col_values(10)
    headers.append(cvir_plus.pop(0))
    cvir_minus = sheet1.col_values(11)
    headers.append(cvir_minus.pop(0))
    mvir = sheet1.col_values(12)
    headers.append(mvir.pop(0))
    mvir_plus = sheet1.col_values(13)
    headers.append(mvir_plus.pop(0))
    mvir_minus = sheet1.col_values(14)
    headers.append(mvir_minus.pop(0))
    short_refs = sheet1.col_values(15)
    headers.append(short_refs.pop(0))
    orig_convention = sheet1.col_values(16)
    headers.append(orig_convention.pop(0))
    cosmology = sheet1.col_values(17)
    headers.append(cosmology.pop(0))
    return clusters,redshift,methods,c200,c200_plus,c200_minus,m200,m200_plus,m200_minus,cvir,cvir_plus,cvir_minus,mvir,mvir_plus,mvir_minus,short_refs,orig_convention,cosmology

def startup_old():
    # Opening Excel File
    wb = open_workbook('/Users/groenera/Desktop/Dropbox/Private/Research/DataFiles/ObservedClusterConcsDB/cm_data.xlsx')

    # Creating Data Structures
    clusters,redshifts,methods,refs = ([],[],[],[])
    origdelta,cosmology,notes = ([],[],[])
    c200,c200_p,c200_m = ([],[],[])
    m200,m200_p,m200_m = ([],[],[])
    cvir,cvir_p,cvir_m = ([],[],[])
    mvir,mvir_p,mvir_m = ([],[],[])

    # Filling Data Structures From Worksheet
    for sheet in wb.sheets():
        if sheet.name == 'Sheet1':
            for i in range(sheet.nrows):
                if i != 0:
                    clusters.append("{}".format(sheet.row(i)[0]).split("'")[1])
                    redshifts.append(np.float("{}".format(sheet.row(i)[1]).split(':')[1]))
                    methods.append("{}".format(sheet.row(i)[2]).split("'")[1])
                    refs.append("{}".format(sheet.row(i)[15]).split("'")[1])
                    origdelta.append("{}".format(sheet.row(i)[-3]).split(":")[1])
                    cosmology.append("{}".format(sheet.row(i)[17]).split("'")[1])
                    notes.append("{}".format(sheet.row(i)[-1]).split("'")[1])
                    # c200
                    try:
                        c200.append(np.float("{}".format(sheet.row(i)[3]).split(':')[1]))
                    except ValueError:
                        c200.append("{}".format(sheet.row(i)[3]).split(':')[1])
                    try:
                        c200_p.append(np.float("{}".format(sheet.row(i)[4]).split(':')[1]))
                    except ValueError:
                        c200_p.append("{}".format(sheet.row(i)[4]).split(':')[1])
                    try:
                        c200_m.append(np.float("{}".format(sheet.row(i)[5]).split(':')[1]))
                    except ValueError:
                        c200_m.append("{}".format(sheet.row(i)[5]).split(':')[1])
                    # m200
                    try:
                        m200.append(np.float("{}".format(sheet.row(i)[6]).split(':')[1]))
                    except ValueError:
                        m200.append("{}".format(sheet.row(i)[6]).split(':')[1])
                    try:
                        m200_p.append(np.float("{}".format(sheet.row(i)[7]).split(':')[1]))
                    except ValueError:
                        m200_p.append("{}".format(sheet.row(i)[7]).split(':')[1])
                    try:
                        m200_m.append(np.float("{}".format(sheet.row(i)[8]).split(':')[1]))
                    except ValueError:
                        m200_m.append("{}".format(sheet.row(i)[8]).split(':')[1])
                    # cvir
                    try:
                        cvir.append(np.float("{}".format(sheet.row(i)[9]).split(':')[1]))
                    except ValueError:
                        cvir.append("{}".format(sheet.row(i)[9]).split(':')[1])
                    try:
                        cvir_p.append(np.float("{}".format(sheet.row(i)[10]).split(':')[1]))
                    except ValueError:
                        cvir_p.append("{}".format(sheet.row(i)[10]).split(':')[1])
                    try:
                        cvir_m.append(np.float("{}".format(sheet.row(i)[11]).split(':')[1]))
                    except ValueError:
                        cvir_m.append("{}".format(sheet.row(i)[11]).split(':')[1])
                    # mvir
                    try:
                        mvir.append(np.float("{}".format(sheet.row(i)[12]).split(':')[1]))
                    except ValueError:
                        mvir.append("{}".format(sheet.row(i)[12]).split(':')[1])
                    try:
                        mvir_p.append(np.float("{}".format(sheet.row(i)[13]).split(':')[1]))
                    except ValueError:
                        mvir_p.append("{}".format(sheet.row(i)[13]).split(':')[1])
                    try:
                        mvir_m.append(np.float("{}".format(sheet.row(i)[14]).split(':')[1]))
                    except ValueError:
                        mvir_m.append("{}".format(sheet.row(i)[14]).split(':')[1])
                
                
    # Converting nans to np.nans, and dealing with infinities
    for i in range(len(c200)):
        # c200
        if type(c200[i]) == str:
            if c200[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif c200[i].split("'")[1] == 'nan':
                c200[i] = np.nan
        if type(c200_p[i]) == str:
            if c200_p[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif c200_p[i].split("'")[1] == 'nan':
                c200_p[i] = np.nan
            elif c200_p[i].split("'")[1] == 'infty':
                c200_p[i] = 1e6
                print "Infty in c200_p: {}".format(clusters[i])
        if type(c200_m[i]) == str:
            if c200_m[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif c200_m[i].split("'")[1] == 'nan':
                c200_m[i] = np.nan
            elif c200_m[i].split("'")[1] == 'infty':
                c200_m[i] = 1e6
                print "Infty in c200_m: {}".format(clusters[i])
        # m200
        if type(m200[i]) == str:
            if m200[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif m200[i].split("'")[1] == 'nan':
                m200[i] = np.nan
        if type(m200_p[i]) == str:
            if m200_p[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif m200_p[i].split("'")[1] == 'nan':
                m200_p[i] = np.nan
            elif m200_p[i].split("'")[1] == 'infty':
                m200_p[i] = 1e6
                print "Infty in m200_p: {}".format(clusters[i])     
        if type(m200_m[i]) == str:
            if m200_m[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif m200_m[i].split("'")[1] == 'nan':
                m200_m[i] = np.nan
            elif m200_m[i].split("'")[1] == 'infty':
                m200_m[i] = 1e6
                print "Infty in m200_m: {}".format(clusters[i])     
        # cvir
        if type(cvir[i]) == str:
            if cvir[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif cvir[i].split("'")[1] == 'nan':
                cvir[i] = np.nan
        if type(cvir_p[i]) == str:
            if cvir_p[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif cvir_p[i].split("'")[1] == 'nan':
                cvir_p[i] = np.nan
            elif cvir_p[i].split("'")[1] == 'infty':
                cvir_p[i] = 1e6
                print "Infty in cvir_p: {}".format(clusters[i])     
        if type(cvir_m[i]) == str:
            if cvir_m[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif cvir_m[i].split("'")[1] == 'nan':
                cvir_m[i] = np.nan
            elif cvir_m[i].split("'")[1] == 'infty':
                cvir_m[i] = 1e6
                print "Infty in cvir_m: {}".format(clusters[i])
        # mvir
        if type(mvir[i]) == str:
            if mvir[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif mvir[i].split("'")[1] == 'nan':
                mvir[i] = np.nan
        if type(mvir_p[i]) == str:
            if mvir_p[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif mvir_p[i].split("'")[1] == 'nan':
                mvir_p[i] = np.nan
            elif mvir_p[i].split("'")[1] == 'infty':
                mvir_p[i] = 1e6
                print "Infty mvir_p: {}".format(clusters[i])
        if type(mvir_m[i]) == str:
            if mvir_m[i].split("'")[1] == 'TBD':
                #print "To be determined"
                continue
            elif mvir_m[i].split("'")[1] == 'nan':
                mvir_m[i] = np.nan
            elif mvir_m[i].split("'")[1] == 'infty':
                mvir_m[i] = 1e6
                print "Infty in mvir_m: {}".format(clusters[i])
                
def scatter_flag(index,tmp_m,tmp_m_p,tmp_m_m,tmp_c,tmp_c_p,tmp_c_m):
    # nans may be passed in a string. convert to np.nan.
    if tmp_m_p[index] == u'nan':
        tmp_m_p[index] = np.nan
    if tmp_m_m[index] == u'nan':
        tmp_m_m[index] = np.nan
    if tmp_c_p[index] == u'nan':
        tmp_c_p[index] = np.nan
    if tmp_c_m[index] == u'nan':
        tmp_c_m[index] = np.nan
    
    # check first that m and c are present and floating point numbers
    if type(tmp_m[index])!=str and type(tmp_c[index]):
        # 1 means it's a np.nan, 0 means it's not
        check = [np.isnan(tmp_m_p[index]),np.isnan(tmp_m_m[index]),
                 np.isnan(tmp_c_p[index]),np.isnan(tmp_c_m[index])]
        # check if all there
        if check == [0,0,0,0]:
            return 0
        # check if none are there
        if check == [1,1,1,1]:
            return 1
        # only mass errorbars
        if check == [0,0,1,1]:
            return 2
        # only conc errorbars
        if check == [1,1,0,0]:
            return 3
        # mass upper & conc upper
        if check == [0,1,0,1]:
            return 4
        # mass lower & conc lower
        if check == [1,0,1,0]:
            return 5
        # mass upper & conc lower
        if check == [0,1,1,0]:
            return 6
        # mass lower & conc upper
        if check == [1,0,0,1]:
            return 7
    else:
        print "Bad conc or mass value."

def plotwithflag(ax,flag,col,tmp_m,tmp_m_p,tmp_m_m,tmp_c,tmp_c_p,tmp_c_m,tmp_z,witherrors):
    if flag == 0:
        if witherrors == True:
            #ipdb.set_trace()
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([abs(tmp_c_m)],[tmp_c_p]),
                        xerr=([abs(tmp_m_m)],[tmp_m_p]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')            
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
    if flag == 1:
        ax.scatter(tmp_m,
                    (1+tmp_z)*tmp_c,
                    color=col,marker='.')
    if flag == 2:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        xerr=([abs(tmp_m_m)],[tmp_m_p]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
    if flag == 3:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([abs(tmp_c_m)],[tmp_c_p]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
    if flag == 4:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([0],[tmp_c_p]),
                        xerr=([0],[tmp_m_p]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
    if flag == 5:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([tmp_c_m],[0]),
                        xerr=([tmp_m_m],[0]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
    if flag == 6:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([tmp_c_m],[0]),
                        xerr=([0],[tmp_m_p]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(np.log10(tmp_m),
                        np.log10((1+tmp_z)*tmp_c),
                        color=col,marker='.')
    if flag == 7:
        if witherrors == True:
            ax.errorbar(tmp_m,
                        (1+tmp_z)*tmp_c,
                        yerr=([0],[tmp_c_p]),
                        xerr=([tmp_m_m],[0]),
                        color=col,
                        elinewidth=2)
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        else:
            ax.scatter(tmp_m,
                        (1+tmp_z)*tmp_c,
                        color=col,marker='.')
        
def colorselect(methods):
    if methods == 'X-ray':
        col = 'green'
    elif methods == 'WL':
        col='purple'
    elif methods == 'SL':
        col = 'red'
    elif methods == 'WL+SL':
        col = 'black'
    elif methods == 'CM':
        col = 'blue'
    elif methods == 'LOSVD':
        col = 'orange'
    return col
        
def scatter_full_sample(delta='200',witherrors=True,coloredmethods=True,method='all',legend=True):
    fig = plt.figure(figsize=(8,8))
    if delta =='200':
        tmp_z,tmp_methods = ([],[])
        tmp_m200,tmp_m200_p,tmp_m200_m = ([],[],[])
        tmp_c200,tmp_c200_p,tmp_c200_m = ([],[],[])
        # the same strings/np.nans should occur in c200 list (if not, something is wrong;
        # I cannot include data unless m200/c200 pair is reported)
        for i in range(len(m200)): 
            if type(m200[i]) != str: # no strings
                if np.isnan(m200[i]) == False: # no np.nans
                    tmp_m200.append(m200[i])
                    tmp_m200_p.append(m200_p[i])
                    tmp_m200_m.append(m200_m[i])
                    tmp_c200.append(c200[i])
                    tmp_c200_p.append(c200_p[i])
                    tmp_c200_m.append(c200_m[i])
                    tmp_z.append(redshifts[i])
                    tmp_methods.append(methods[i])
        ax = plt.subplot(111)
        ax.set_xscale("log")
        ax.set_yscale("log")
        for i in range(len(tmp_m200)):
            if coloredmethods == True:
                col = colorselect(tmp_methods[i])
            else:
                col= 'blue'
            if method in ['all','All']:
                flag = scatter_flag(i,tmp_m200,tmp_m200_p,tmp_m200_m,
                                    tmp_c200,tmp_c200_p,tmp_c200_m)
                plotwithflag(ax,flag,col,witherrors,
                            tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                            tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                            tmp_z[i])
            else:
                if tmp_methods[i].upper() == method.upper():
                    flag = scatter_flag(i,tmp_m200,tmp_m200_p,tmp_m200_m,
                                    tmp_c200,tmp_c200_p,tmp_c200_m)
                    plotwithflag(ax,flag,col,witherrors,
                                tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                                tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                                tmp_z[i])
        # making legend
        if legend is True:
            plt.scatter(1e6,1e6,color='green',label='X-ray')
            plt.scatter(1e6,1e6,color='purple',label='WL')
            plt.scatter(1e6,1e6,color='red',label='SL')
            plt.scatter(1e6,1e6,color='black',label='WL+SL')
            plt.scatter(1e6,1e6,color='blue',label='CM')
            plt.scatter(1e6,1e6,color='orange',label='LOSVD')
            plt.legend(loc=0,scatterpoints=1)
        plt.xlim(0.1,120)
        plt.ylim(0.2,160)
        plt.xlabel(r'$\mathrm{M_{200}}/{10^{14}\mathrm{M} \odot}$')
        plt.ylabel(r'$(1+z)c_{200}$')
        plt.show()
    elif delta in ['vir','virial','Vir','Virial']:
        frame = sys._getframe(1)
        tmp_z,tmp_methods = ([],[])
        tmp_mvir,tmp_mvir_p,tmp_mvir_m = ([],[],[])
        tmp_cvir,tmp_cvir_p,tmp_cvir_m = ([],[],[])
        # the same strings/np.nans should occur in cvir list (if not, something is wrong;
        # I cannot include data unless mvir/cvir pair is reported)
        for i in range(len(mvir)): 
            if type(mvir[i]) != str: # no strings
                if np.isnan(mvir[i]) == False: # no np.nans
                    tmp_mvir.append(mvir[i])
                    tmp_mvir_p.append(mvir_p[i])
                    tmp_mvir_m.append(mvir_m[i])
                    tmp_cvir.append(cvir[i])
                    tmp_cvir_p.append(cvir_p[i])
                    tmp_cvir_m.append(cvir_m[i])
                    tmp_z.append(redshifts[i])
                    tmp_methods.append(methods[i])
        if frame.f_code.co_name != 'scatter_grid':
            ax = plt.subplot(111)
            ax.set_xscale("log")
            ax.set_yscale("log")
        else:
            return tmp_mvir,tmp_mvir_p,tmp_mvir_m,tmp_cvir,tmp_cvir_p,tmp_cvir_m,tmp_z,tmp_methods
        for i in range(len(tmp_mvir)):
            if coloredmethods == True:
                col = colorselect(tmp_methods[i])
            else:
                col= 'blue'
            if method in ['all','All']:
                flag = scatter_flag(i,tmp_mvir,tmp_mvir_p,tmp_mvir_m,
                                    tmp_cvir,tmp_cvir_p,tmp_cvir_m)
                plotwithflag(ax,flag,col,witherrors,
                            tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                            tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                            tmp_z[i])
            else:
                if tmp_methods[i].upper() == method.upper():
                    flag = scatter_flag(i,tmp_mvir,tmp_mvir_p,tmp_mvir_m,
                                    tmp_cvir,tmp_cvir_p,tmp_cvir_m)
                    plotwithflag(ax,flag,col,witherrors,
                                tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                                tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                                tmp_z[i])
        # making legend
        if legend is True:
            plt.scatter(1e6,1e6,color='green',label='X-ray')
            plt.scatter(1e6,1e6,color='purple',label='WL')
            plt.scatter(1e6,1e6,color='red',label='SL')
            plt.scatter(1e6,1e6,color='black',label='WL+SL')
            plt.scatter(1e6,1e6,color='blue',label='CM')
            plt.scatter(1e6,1e6,color='orange',label='LOSVD')
            plt.legend(loc=0,scatterpoints=1)
        plt.xlim(0.1,110)
        plt.ylim(0.1,160)
        plt.xlabel(r'$\mathrm{M_{vir}}/{10^{14}\mathrm{M} \odot}$')
        plt.ylabel(r'$(1+z)c_{vir}$')
        plt.show()

# Kind of deprecated at the moment.
def explore_data():
    # First make sure there are the same number of datapoints for cvir/mvir
    assert len(cvir) == len(mvir)
    # Find list indices for cvir/mvir if type is string
    cvir_strings = [i for i in range(len(cvir)) if type(cvir[i]) == str]
    mvir_strings = [i for i in range(len(mvir)) if type(mvir[i]) == str]
    # Make sure the list elements (for strings) for cvir/mvir are the same
    assert cvir_strings == mvir_strings
    # Keep only non-strings for cvir/mvir
    cvir_vals = [cvir[i] for i in range(len(cvir)) if i not in cvir_strings]
    mvir_vals = [mvir[i] for i in range(len(mvir)) if i not in mvir_strings]
    # Make sure that there are references for every cvir measurement
    assert len(cvir) == len(refs)
    # Keep refs for only float values of cvir/mvir
    refs_nostr = [refs[i] for i in range(len(refs)) if i not in cvir_strings]
    # Make sure that there are the same number of datapoints for cvir_p/cvir_m, mvir_p/mvir_m
    assert len(cvir_p) == len(cvir_m)
    assert len(cvir_p) == len(cvir)
    assert len(mvir_p) == len(mvir_m)
    assert len(mvir_p) == len(mvir)
    # Find list indices for cvir_p/cvir_m, mvir_p/mvir_m if type is string
    cvir_p_strings = [i for i in range(len(cvir_p)) if type(cvir_p[i]) == str]
    cvir_m_strings = [i for i in range(len(cvir_m)) if type(cvir_m[i]) == str]
    mvir_p_strings = [i for i in range(len(mvir_p)) if type(mvir_p[i]) == str]
    mvir_m_strings = [i for i in range(len(mvir_m)) if type(mvir_m[i]) == str]
    # Make sure the list elements (for strings) for cvir_p/cvir_m, mvir_p/mvir_m are the same
    assert cvir_p_strings == cvir_m_strings
    assert cvir_p_strings == cvir_strings
    assert mvir_p_strings == mvir_m_strings
    assert mvir_p_strings == mvir_strings
    # Keep only non-strings for for cvir_p/cvir_m, mvir_p/mvir_m
    cvir_p_vals = [cvir_p[i] for i in range(len(cvir)) if i not in cvir_strings]
    cvir_m_vals = [cvir_m[i] for i in range(len(cvir)) if i not in cvir_strings]
    mvir_p_vals = [mvir_p[i] for i in range(len(cvir)) if i not in mvir_strings]
    mvir_m_vals = [mvir_m[i] for i in range(len(cvir)) if i not in mvir_strings]    
    # Find list indices for cvir_p/cvir_m, mvir_p/mvir_m if type is np.nan
    cvir_p_nans = [i for i in range(len(cvir_p_vals)) if np.isnan(cvir_p_vals[i])]
    cvir_m_nans = [i for i in range(len(cvir_m_vals)) if np.isnan(cvir_m_vals[i])]
    mvir_p_nans = [i for i in range(len(mvir_p_vals)) if np.isnan(mvir_p_vals[i])]
    mvir_m_nans = [i for i in range(len(mvir_m_vals)) if np.isnan(mvir_m_vals[i])]
    # Make sure same amount of nans for each uncertainty
    assert cvir_p_nans == cvir_m_nans
    assert mvir_p_nans == mvir_m_nans
    # Make sure only floats are left
    cvir_p_vals = [cvir_p_vals[i] for i in range(len(cvir_p_vals)) if i not in cvir_p_nans]
    cvir_m_vals = [cvir_m_vals[i] for i in range(len(cvir_m_vals)) if i not in cvir_m_nans]
    mvir_p_vals = [mvir_p_vals[i] for i in range(len(mvir_p_vals)) if i not in mvir_p_nans]
    mvir_m_vals = [mvir_m_vals[i] for i in range(len(mvir_m_vals)) if i not in mvir_m_nans]
    ipdb.set_trace()
    # get all refs which are now associated with    
    refs_nonan = [refs_nostr[i] for i in range(len(refs_nostr)) if i not in cvir_p_nans]
    uncertainty_diffs = [abs(abs(cvir_p_vals[i])-abs(cvir_m_vals[i])) for i in range(len(cvir_p_vals))]
    print "Out of a total of {} halos with concentration uncertainties reported, {} of them have symmetric errorbars.".format(len(cvir_p_vals),len([i for i in range(len(cvir_p_vals)) if uncertainty_diffs[i] == 0]))
    asym_refs = [refs_nonan[i] for i in range(len(refs_nonan)) if uncertainty_diffs[i] != 0]
    print "References which report assymetric concentration uncertainties are:"
    reftable = Table()
    unique_refs = [i for i in iter(set(asym_refs))]
    reftable['Ref.'] = unique_refs
    reftable['# Asymmetric Clusters'] = [refs_nonan.count(i) for i in unique_refs]
    reftable['# of Total Clusters in Paper'] = [refs.count(i) for i in unique_refs]
    reftable.pprint()
    ipdb.set_trace()

def summarize_refs():
    clusters,redshift,methods,c200,c200_p,c200_m,m200,m200_p,m200_m,cvir,cvir_p,cvir_m,mvir,mvir_p,mvir_m,refs,orig_convention,cosmology = startup()
    unique_refs = [i for i in iter(set(refs))]
    unique_refs_count = [refs.count(i) for i in unique_refs]
    plt.hist(unique_refs_count,bins=range(max(unique_refs_count)+1))
    plt.xlabel(r'$\mathrm{N_{cl}}$/Paper',fontsize=18)
    plt.ylabel(r'$\mathrm{N_{papers}}$',fontsize=18)
    largest = unique_refs[unique_refs_count.index(max(unique_refs_count))]
    sortedrefscount = copy.deepcopy(unique_refs_count)
    sortedrefscount.sort()
    #ipdb.set_trace()
    secondlargest = unique_refs[unique_refs_count.index(sortedrefscount[-2])] 
    thirdlargest = unique_refs[unique_refs_count.index(sortedrefscount[-3])]
    plt.text(max(unique_refs_count)-3,2,'{}'.format(largest))
    plt.text(sortedrefscount[-2]-3,2,'{}'.format(secondlargest))
    plt.text(sortedrefscount[-3]-3,2,'{}'.format(thirdlargest))
    plt.show()

def summarize_methods():
    clusters,redshift,methods,c200,c200_p,c200_m,m200,m200_p,m200_m,cvir,cvir_p,cvir_m,mvir,mvir_p,mvir_m,refs,orig_convention,cosmology = startup()
    unique_methods = [i for i in iter(set(methods))]
    unique_method_count = [methods.count(i) for i in unique_methods]
    t = Table()
    t['Method'] = unique_methods
    t['N. Meas.'] = unique_method_count
    ascii.write(t,format='latex')

def popular_clusters():
    clusters,redshift,methods,c200,c200_p,c200_m,m200,m200_p,m200_m,cvir,cvir_p,cvir_m,mvir,mvir_p,mvir_m,refs,orig_convention,cosmology = startup()
    unique_clusters = [i for i in iter(set(clusters))]
    unique_cluster_count = [clusters.count(i) for i in unique_clusters]
    multiple_cl = [unique_clusters[i] for i in range(len(unique_clusters)) if unique_cluster_count[i] > 1]
    multiple_ct = [unique_cluster_count[i] for i in range(len(unique_clusters)) if unique_cluster_count[i] > 1]
    multiple_methods = [[i for i in iter(set([methods[i] for i, x in enumerate(clusters) if x == j]))] for j in multiple_cl]
    multiple_refs = [[i for i in iter(set([refs[i] for i, x in enumerate(clusters) if x == j]))] for j in multiple_cl]
    t = Table()
    t['Cluster'] = multiple_cl
    t['N. Meas.'] = multiple_ct
    t['Method(s)'] = multiple_methods
    t['Refs.'] = multiple_refs
    ascii.write(t,format='latex')

def uncertainty_summary(redshift_plot = False):
    clusters,redshift,methods,c200,c200_p,c200_m,m200,m200_p,m200_m,cvir,cvir_p,cvir_m,mvir,mvir_p,mvir_m,refs,orig_convention,cosmology = startup()
    # c/m measurements first
    # first figure out how many have c/M measurements; gets rid of the TBD entries
    complete_cvir = [i for i in range(len(clusters)) if type(cvir[i]) != str]
    complete_mvir = [i for i in range(len(clusters)) if type(mvir[i]) != str]
    assert complete_cvir == complete_mvir # same TBD entries for c and M
    # filter out nans
    complete_cvir = [i for i in complete_cvir if cvir[i] not in [u'nan']]
    complete_mvir = [i for i in complete_mvir if mvir[i] not in [u'nan']]
    match = []
    # pop elements off of complete_cvir if also in complete_mvir; destroys complete_cvir
    for i in complete_mvir:
        match.append(complete_cvir.pop(complete_cvir.index(i)))
    missing_mass = complete_cvir
    complete_cvir = match

    # find papers which only report concentration
    tmp_refs = [refs[i] for i in missing_mass]
    missing_mass_refs = [i for i in iter(set(tmp_refs))]
    missing_mass_refs_ct = [tmp_refs.count(i) for i in missing_mass_refs]
    total_refs_ct = [refs.count(i) for i in missing_mass_refs]
    t = Table()
    t['Refs.'] = missing_mass_refs
    t['No. Missing Mass Clusters'] = missing_mass_refs_ct
    t['Tot. No. Clusters in Ref.'] = total_refs_ct
    ascii.write(t,format='latex')

    # of the c/M measurement pairs, how many have both errorbars?
    complete_cvir_p = [cvir_p[i] for i in match]
    complete_cvir_m = [cvir_m[i] for i in match]
    complete_mvir_p = [mvir_p[i] for i in match]
    complete_mvir_m = [mvir_m[i] for i in match]
    allfour = [i for i in range(len(match)) if
               complete_cvir_p[i] not in [u'nan'] and
               complete_cvir_m[i] not in [u'nan'] and
               complete_mvir_p[i] not in [u'nan'] and
               complete_mvir_m[i] not in [u'nan']]
    none = [i for i in range(len(match)) if
                complete_cvir_p[i] in [u'nan'] and
                complete_cvir_m[i] in [u'nan'] and
                complete_mvir_p[i] in [u'nan'] and
                complete_mvir_m[i] in [u'nan']]
    justconc = [i for i in range(len(match)) if
                complete_cvir_p[i] not in [u'nan'] and
                complete_cvir_m[i] not in [u'nan'] and
                complete_mvir_p[i] in [u'nan'] and
                complete_mvir_m[i] in [u'nan']]
    justmass = [i for i in range(len(match)) if
                complete_cvir_p[i] in [u'nan'] and
                complete_cvir_m[i] in [u'nan'] and
                complete_mvir_p[i] not in [u'nan'] and
                complete_mvir_m[i] not in [u'nan']]
    allfour_sym = [i for i in allfour
                   if abs(complete_cvir_p[i]) == abs(complete_cvir_m[i])
                   and abs(complete_mvir_p[i]) == abs(complete_mvir_m[i])]
    allfour_asym = [i for i in allfour if i not in allfour_sym]
    justconc_sym = [i for i in justconc
                   if abs(complete_cvir_p[i]) == abs(complete_cvir_m[i])]
    justconc_asym = [i for i in justconc if i not in justconc_sym]
    justmass_sym = [i for i in justmass
                   if abs(complete_mvir_p[i]) == abs(complete_mvir_m[i])]
    justmass_asym = [i for i in justmass if i not in justmass_sym]
    t2 = Table()
    t2['Uncertainty Type'] = ['Both c/M','Just M','Just c','Neither']
    t2['Symmetric'] = [len(allfour_sym),len(justmass_sym),len(justconc_sym),'--']
    t2['Asymmetric'] = [len(allfour_asym),len(justmass_asym),len(justconc_asym),'--']
    t2['Total'] = [len(allfour),len(justmass),len(justconc),len(none)]
    ascii.write(t2,format='latex')

    # redshift distribution for clusters with both c/m measurements
    if redshift_plot == True:
        match_redshifts = [redshifts[i] for i in match]
        plt.hist(match_redshifts, histtype='stepfilled')
        plt.xlabel('redshift')
        plt.ylabel('Number')
        plt.show()

def scatter_uncertainty_sample(delta='200',witherrors=True,coloredmethods=True,method='all'):
    fig = plt.figure(figsize=(8,8))
    if delta =='200':
        tmp_z,tmp_methods = ([],[])
        tmp_m200,tmp_m200_p,tmp_m200_m = ([],[],[])
        tmp_c200,tmp_c200_p,tmp_c200_m = ([],[],[])
        # the same strings/np.nans should occur in c200 list (if not, something is wrong;
        # I cannot include data unless m200/c200 pair is reported)
        for i in range(len(m200)): 
            if type(m200[i]) != str: # no strings
                if np.isnan(m200[i]) == False: # no np.nans
                    tmp_m200.append(m200[i])
                    tmp_m200_p.append(m200_p[i])
                    tmp_m200_m.append(m200_m[i])
                    tmp_c200.append(c200[i])
                    tmp_c200_p.append(c200_p[i])
                    tmp_c200_m.append(c200_m[i])
                    tmp_z.append(redshifts[i])
                    tmp_methods.append(methods[i])
        allfour = 0
        none = 0
        ax = plt.subplot(111)
        ax.set_xscale("log")
        ax.set_yscale("log")
        for i in range(len(tmp_m200)):
            if coloredmethods == True:
                col = colorselect(tmp_methods[i])
            else:
                col= 'blue'
            if method in ['all','All']:
                flag = scatter_flag(i,tmp_m200,tmp_m200_p,tmp_m200_m,
                                    tmp_c200,tmp_c200_p,tmp_c200_m)
                if flag in [0,1]:
                    if flag == 0:
                        if abs(tmp_m200_p[i]) == abs(tmp_m200_m[i]) and abs(tmp_c200_p[i]) == abs(tmp_c200_m[i]):
                            plotwithflag(ax,flag,col,witherrors,
                                        tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                                        tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                                        tmp_z[i])
                    if flag == 1:
                        plotwithflag(ax,flag,col,witherrors,
                                        tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                                        tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                                        tmp_z[i])
            else:
                if tmp_methods[i].upper() == method.upper():
                    flag = scatter_flag(i,tmp_m200,tmp_m200_p,tmp_m200_m,
                                    tmp_c200,tmp_c200_p,tmp_c200_m)
                    if flag in [0,1]:
                        if flag == 0:
                            if abs(tmp_m200_p[i]) == abs(tmp_m200_m[i]) and abs(tmp_c200_p[i]) == abs(tmp_c200_m[i]):
                                plotwithflag(ax,flag,col,witherrors,
                                            tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                                            tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                                            tmp_z[i])
                        if flag == 1:
                            plotwithflag(ax,flag,col,witherrors,
                                        tmp_m200[i],tmp_m200_p[i],tmp_m200_m[i],
                                        tmp_c200[i],tmp_c200_p[i],tmp_c200_m[i],
                                        tmp_z[i])
        plt.xlim(0.1,110)
        plt.ylim(0.1,160)
        plt.xlabel(r'$\mathrm{M_{200}}/{10^{14}\mathrm{M} \odot}$')
        plt.ylabel(r'$(1+z)c_{200}$')
        plt.show()
    elif delta in ['vir','virial','Vir','Virial']:
        tmp_z,tmp_methods = ([],[])
        tmp_mvir,tmp_mvir_p,tmp_mvir_m = ([],[],[])
        tmp_cvir,tmp_cvir_p,tmp_cvir_m = ([],[],[])
        # the same strings/np.nans should occur in cvir list (if not, something is wrong;
        # I cannot include data unless mvir/cvir pair is reported)
        for i in range(len(mvir)): 
            if type(mvir[i]) != str: # no strings
                if np.isnan(mvir[i]) == False: # no np.nans
                    tmp_mvir.append(mvir[i])
                    tmp_mvir_p.append(mvir_p[i])
                    tmp_mvir_m.append(mvir_m[i])
                    tmp_cvir.append(cvir[i])
                    tmp_cvir_p.append(cvir_p[i])
                    tmp_cvir_m.append(cvir_m[i])
                    tmp_z.append(redshifts[i])
                    tmp_methods.append(methods[i])
        ax = plt.subplot(111)
        ax.set_xscale("log")
        ax.set_yscale("log")
        allfour = 0
        none = 0
        for i in range(len(tmp_mvir)):
            if coloredmethods == True:
                col = colorselect(tmp_methods[i])
            else:
                col= 'blue'
            if method in ['all','All']:
                flag = scatter_flag(i,tmp_mvir,tmp_mvir_p,tmp_mvir_m,
                                    tmp_cvir,tmp_cvir_p,tmp_cvir_m)
                if flag in [0,1]:
                    if flag == 0:
                        if abs(tmp_mvir_p[i]) == abs(tmp_mvir_m[i]) and abs(tmp_cvir_p[i]) == abs(tmp_cvir_m[i]):
                            plotwithflag(ax,flag,col,witherrors,
                                        tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                                        tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                                        tmp_z[i])
                            allfour += 1
                    if flag == 1:
                        plotwithflag(ax,flag,col,witherrors,
                                        tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                                        tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                                        tmp_z[i])
                        none += 1
            else:
                if tmp_methods[i].upper() == method.upper():
                    flag = scatter_flag(i,tmp_mvir,tmp_mvir_p,tmp_mvir_m,
                                    tmp_cvir,tmp_cvir_p,tmp_cvir_m)
                    if flag in [0,1]:
                        if flag == 0:
                            if abs(tmp_m200_p[i]) == abs(tmp_m200_m[i]) and abs(tmp_c200_p[i]) == abs(tmp_c200_m[i]):
                                plotwithflag(ax,flag,col,witherrors,
                                            tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                                            tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                                            tmp_z[i])
                        if flag == 1:
                            plotwithflag(ax,flag,col,witherrors,
                                        tmp_mvir[i],tmp_mvir_p[i],tmp_mvir_m[i],
                                        tmp_cvir[i],tmp_cvir_p[i],tmp_cvir_m[i],
                                        tmp_z[i])
        plt.xlim(0.1,110)
        plt.ylim(0.1,160)
        plt.xlabel(r'$\mathrm{M_{vir}}/{10^{14}\mathrm{M} \odot}$')
        plt.ylabel(r'$(1+z)c_{vir}$')
        plt.show()

def check_cosmology(plot_cosmo=False):
    clusters,redshift,methods,c200,c200_p,c200_m,m200,m200_p,m200_m,cvir,cvir_p,cvir_m,mvir,mvir_p,mvir_m,refs,orig_convention,cosmology = startup()
    method_list = []
    cosmo_list = []
    for i in range(len(cosmology)):
        tmp = cosmology[i].strip('(').strip(')').split('/')
        if tmp[0] != '0.3' or tmp[1] != '0.7':
            method_list.append(methods[i])
            if cosmology[i] != 'TBD':
                cosmo_list.append(cosmology[i])
        else:
            print cosmology[i]
    if plot_cosmo is True:
        omegam = [i.strip('(').strip(')').split('/')[0] for i in cosmo_list]
        omegal = [i.strip('(').strip(')').split('/')[1] for i in cosmo_list]
        for i in range(len(cosmo_list)):
            plt.scatter(omegam[i],omegal[i])
        plt.xlabel(r'$\Omega_{m}$')
        plt.ylabel(r'$\Omega_{\Lambda}$')
        plt.xlim(0,1)
        plt.ylim(0,1)
        plt.show()
    return method_list,cosmo_list

def cmrelation_co07_1(mvir_min,mvir_max,z,c0=14.8,alpha=-0.14,mstar=1.3e13/0.7):
    mlist = np.linspace(mvir_min,mvir_max,1000)
    clist = [(c0/(1+z)) * ((mlist[i]/mstar)**(alpha)) for i in range(len(mlist))]
    return mlist, clist, z

def cmrelation_bu01_1(mvir_min,mvir_max,z,c0=9.0,alpha=-0.13,mstar=1.3e13/0.7):
    mlist = np.linspace(mvir_min,mvir_max,1000)
    clist = [(c0/(1+z)) * ((mlist[i]/mstar)**(alpha)) for i in range(len(mlist))]
    return mlist, clist, z

def cmrelation_he07_1(mvir_min,mvir_max,z,c0=12.3,alpha=-0.13,mstar=1.3e13/0.7):
    mlist = np.linspace(mvir_min,mvir_max,1000)
    clist = [(c0/(1+z)) * ((mlist[i]/mstar)**(alpha)) for i in range(len(mlist))]
    return mlist, clist, z

def cmrelation_gr14_1(mvir_min,mvir_max,z,c0=4.775,alpha=-0.056,mstar=1e14/0.7):
    mlist = np.linspace(mvir_min,mvir_max,1000)
    clist = [(c0/(1+z)) * ((mlist[i]/mstar)**(alpha)) for i in range(len(mlist))]
    return mlist, clist, z

    
def cmrelation_pr11_1(mvir_min,mvir_max,z,omega_m_0=0.3,omega_l_0=0.7):
    def cmin(x):
        return 3.681+(5.033-3.681)*( (1./np.pi)*np.arctan(6.948*(x-0.424)) + 0.5 )
    def siginvmin(x):
        return 1.047 + (1.646-1.047)*( (1./np.pi)*np.arctan(7.386*(x-0.526)) + 0.5 )
    # some preliminary things
    a = 1./(1+z)
    x = (omega_l_0/omega_m_0)**(1./3) * a
    integral = quad(lambda xp: xp**(3./2)/((1+xp**3)**(3./2)),0,x)[0]
    D = (5./2.) * (omega_m_0/omega_l_0)**(1./3) * ((np.sqrt(1+x**3))/(x**(3./2))) * integral
    B0 = cmin(x)/cmin(1.393)
    B1 = siginvmin(x)/siginvmin(1.393)
    # mass dependent quantities
    mlist = np.linspace(mvir_min,mvir_max,1000)
    y = [mlist[i]/1e12/0.7 for i in range(len(mlist))]
    sigma = [D * ((16.9*(y[i]**0.41))/(1+1.102*(y[i]**0.20)+6.22*(y[i]**0.333))) for i in range(len(y))]
    sigmap = [B1 * sigma[i] for i in range(len(sigma))]
    C = [2.881*((sigmap[i]/1.257)**1.022+1)*np.exp(0.060/(sigmap[i]**2)) for i in range(len(sigmap))]
    clist = [B0*C[i] for i in range(len(C))]
    return mlist,clist,z
    
                        
if __name__ == '__main__':
    #scatter_full_sample(delta='vir',witherrors=True,coloredmethods=True,method='all')
    #explore_data()
    #uncertainty_summary(redshift_plot = True)
    #scatter_uncertainty_sample(delta='vir')
    #methods_notLCDM,nonstandard_cosmologies = check_cosmology()
    #ipdb.set_trace()
    '''
    mlistPR11_1_1689,clistPR11_1_1689,zPR11_1_1689 = cmrelation_pr11_1(1e14,2e16,0)
    mlistPR11_1_2137,clistPR11_1_2137,zPR11_1_2137 = cmrelation_pr11_1(1e14,2e16,0.2)
    mlistPR11_1_1835,clistPR11_1_1835,zPR11_1_1835 = cmrelation_pr11_1(1e14,2e16,0.4)
    plt.plot(mlistPR11_1_1689,clistPR11_1_1689)
    plt.plot(mlistPR11_1_2137,clistPR11_1_2137)
    plt.plot(mlistPR11_1_1835,clistPR11_1_1835)
    plt.show()
    '''
    #uncertainty_summary()
    methods, cosmo = check_cosmology(plot_cosmo=True)
    ipdb.set_trace()
