#!/bin/python

import numpy as np
from collections import Counter
from scipy import interpolate
import time
import pymultinest
import warnings
import sys
import argparse

from frb_util_flux import *

dis = AstroDistribution()
cos = Cosmology()
er = EventRate()
tel = Telescope()
lf = Loadfiles()

dnu = 1000.

def lnlik(vpar):
    #global vLOGF, vDME, vLOGW, vN
    #global vSN0, vTs, vG, vBW, vNpol, vFOV, vTime
    #global vLOGF_2d, vDME_2d, vLOGW_2d 
    try:
        norm = np.zeros(vFOV.shape)
        for i in range(len(norm)):
            norm[i] = dis.Norm1D(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i], dnu, vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_fdm = np.zeros(vN.shape)
        for i in range(len(vN)):
            loglik_fdm[i] = np.sum(dis.log_distr_fdmw(dnu, vLOGF_2d[i], vDME_2d[i], vLOGW_2d[i], vpar[1], vpar[2], vpar[3], vpar[4], vpar[5], gtype=fgt)-np.log(norm[i]))
        loglik_norm = np.sum(loglik_fdm)
        rho = np.zeros(vFOV.shape)
        for i in range(len(rho)):
            rho[i] = er.rate_2d(vSN0[i], vBW[i], vNpol[i], vG[i], vTs[i], dnu, vpar[0], vpar[1], vpar[2], vpar[3], vpar[4], vpar[5])
        loglik_poi = np.sum(er.log_dis_poi(rho, vN, vFOV, vTime))
        res = loglik_norm + loglik_poi
        return res
    except:
        print 'Numerical error: @', vpar
        return -1e99

def myprior(cube, ndim, nparams):
    for i in range(ndim):
        cube[i] = vpar_range[i,0]+cube[i] *(vpar_range[i,1]-vpar_range[i,0])
    if bolupper:
        cube[0] = np.power(10., cube[0])
        cube[3] = np.log10(cube[3])

def myloglike(cube, ndim, nparams):
    cube2 = np.zeros(ndim)
    for i in range(0,ndim):
        cube2[i] = cube[i]
    res = lnlik(cube2)
    return res

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='LF Measurements Program for FRB sample')
    parser.add_argument('-fc', action='store', dest='fcat', type=str, help='Input FRB catalog file')
    parser.add_argument('-fs', action='store', dest='fsvy', type=str, help='Input Survey Information file')
    parser.add_argument('-o', action='store', dest='fout', type=str, help='Save the output file')
    parser.add_argument('-g', action='store', dest='fgt', type=str, help='Host galaxy type')
    parser.add_argument('-upper', dest='bolupper', type=bool, help='Bool option: choosing flat prior')
    parser.add_argument('-halo', dest='bolhalo', type=bool, help='Bool option: removing the DM from dark halo')
    
    args = parser.parse_args()
    fcat = args.fcat
    fsvy = args.fsvy
    fout = args.fout
    fgt = args.fgt
    bolupper = args.bolupper
    bolhalo = args.bolhalo
    
    frb_cat = lf.LoadCatalogue(fcat)
    svy_info = lf.LoadSvyInfo(fsvy)
    
    if fgt.find('NE2001')>=0:
        vDME = frb_cat['DM'] - frb_cat['DM_NE2001']
    else:
        vDME = frb_cat['DM'] - frb_cat['DM_YMW16']
    
    if bolhalo:
        vDME = vDME-30.
    
    if fgt.find('ETG')>=0:
        fgt='ETG'
    
    vLOGF = np.log10(frb_cat['S'])
    vLOGW = np.log10(frb_cat['W'])
    vSVY = frb_cat['SURVEY']
    cts = Counter(vSVY)
    
    unique_svys = ['CHIME_1', 'CHIME_2']
    vN = np.array([cts[s] for s in unique_svys])
    
    vSN0 = svy_info['SN0']
    vTs = svy_info['Tsys']
    vG = svy_info['Gain']
    vBW = svy_info['BW']
    vNpol = svy_info['Npol']
    vFOV = svy_info['FOV'] 
    vTime = svy_info['TIME']
    
    vLOGF_2d = [vLOGF[vSVY==s] for s in unique_svys]
    vDME_2d = [vDME[vSVY==s] for s in unique_svys]
    vLOGW_2d = [vLOGW[vSVY==s] for s in unique_svys]
    
    if bolupper:
        # Prioir for Flux Function: 
        # P0: log_alpha (index), P1: log_fs (pseudo flux cut-off), P2: log_f0 (min flux)
        # However, the physical meaning of alpha, logls, logl0 -> alpha, logfs, logf0
        # Let's adjust range to fit flux -2 to 4
        vpara=np.array([-3.0, -3.0, -4.0, 1e2, -1, 0.1])
        vparb=np.array([4.47, 1.1, 4.0, 1e7, 2, 1.0])
    else:
        vpara=np.array([1e-3, -3.0, -4.0, -2.0, -1, 0.1])
        vparb=np.array([3e4, 1.1, 4.0, 4.0, 2, 1.0])
    
    vpar_range=np.dstack((vpara.transpose(),vparb.transpose()))[0,:,:]
    
    
    print '------------par range-----------'
    print vpar_range
    a1=time.clock()
    print myloglike(vpara, len(vpara),len(vpara))
    a2=time.clock()
    print a1,a2
    print "Running Nest Sampling ..."
    # run MultiNest
    pymultinest.run(myloglike, myprior, len(vpara),
                    importance_nested_sampling = False,
                    resume = False,
                    verbose = True,
                    sampling_efficiency = 'model',
                    n_live_points = 1000,
                    outputfiles_basename='nest_out/samp/'+fout)
