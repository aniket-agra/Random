#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 14:57:41 2019

@author: aagrawal
"""

import numpy as np
#from scipy.integrate import quad
#import hankel
from hankel import HankelTransform
#from scipy.interpolate import InterpolatedUnivariateSpline as spline
import fit_profile as fp
import nbodykit.cosmology.correlation as corr
 
popt = fp.main()
keyarr = ('rhoi', 'rhox', 'rs', 'rt', 'be', 'se', 'al', 'gam')
params = dict(zip(keyarr, popt))
epspar = {'rhox' : 1, 'rhoi' : 20, 'rs' : 0.01, 'rt' : 0.01, 'be' : 0.01, \
          'se' : 0.01, 'al' : 0.01, 'gam' : 0.1}
rparr = np.logspace(-1, 1, 20)    
lparams = len(params)
lrparr = len(rparr)
basearr = np.zeros((1, lrparr))
derarr = []
basepar = params.copy()
keyarr = params.keys()
basexi = lambda r : fp.rho_r_parr(r, basepar)


#covm_w_der = np.zeros((lparams+1, lrparr, lrparr))
#pkarr = np.empty([lparams, 1], dtype = 'object')
#ht1 = HankelTransform(nu = 1/2, h = 1e-3)
#ht2 = HankelTransform(nu = -1/2, h = 1e-3)

#basexi = lambda r : np.sqrt(r)*fp.rho_r_parr(r, basepar)
#basepk = lambda kq : np.sqrt(8*np.pi**3/kq)*ht1.transform(basexi, kq, ret_err=False, ret_cumsum=False)

#basepk = lambda kq : np.sqrt(8*np.pi**3/kq)*ht1.transform(basexi, kq, ret_err=False, ret_cumsum=False)
#
#for key in keyarr : 
#    params[key] = basepar[key] + epspar[key]
#    derp = lambda r : fp.rho_r_parr(r, params)
#    params[key] = basepar[key] - epspar[key]
#    derm = lambda r : fp.rho_r_parr(r, params)
#    xirp = lambda r : np.sqrt(r)*(0.5/epspar[key])*(derp(r)-derm(r))
#    pkp = lambda kq : np.sqrt(8*np.pi**3/kq)*ht1.transform(xirp, kq, ret_err=False, ret_cumsum=False)
#    derarr.append(pkp)
#
#kmin = 1e-3 
#karr = np.logspace(-0.7, 1.3, 231)
#for i in np.arange(lparams+1) : 
#    for j in np.arange(lrparr) : 
#        r1 = rparr[j]
#        for l in np.arange(j, lrparr) : 
#            r2 = rparr[l]
#            freqm = np.abs(r1-r2)
#            sqrtm = norm*np.sqrt(freqm)
#            freqp = r1+r2
#            sqrtp = norm*np.sqrt(freqp)
#            if i == 0 : 
#                inte = lambda kq : basepk(kq+kmin)**2/np.sqrt(kq)
#            else : 
#                inte = lambda kq : 2*derarr[i-1](kq+kmin)*basepk(kq+kmin)/np.sqrt(kq)
#            if j == l : 
#                covm_w_der[i, j, l] = quad(basepk, kmin, np.inf, limit = 100)[0] \
#                                    - sqrtp * ht2.transform(inte, freqp, ret_err=False)            
#            else : 
#                covm_w_der[i, j, l] = sqrtm * ht2.transform(inte, freqm, ret_err=False) \
#                                    - sqrtp * ht2.transform(inte, freqp, ret_err=False)
#            covm_w_der[:, j, l] = covm_w_der[:, j, l]/(4*np.pi**2*r1*r2)
#            covm_w_der[:, l, j] = covm_w_der[:, j, l]        



#for i in np.arange(lparams) : 
#   cicon = np.matmul(covi, covm_w_der[i+1, :, :])
#   for j in np.arange(i, lparams) : 
#       cjcon = np.matmul(covi, covm_w_der[j+1, :, :])
#       meander = np.outer(derarr[j, :], derarr[i, :]) + np.outer(derarr[i, :], derarr[j, :])
#       fishmat[i, j] = 0.5*np.trace(np.matmul(cicon, cjcon) + np.matmul(covi, meander))
#       fishmat[j, i] = fishmat[i, j]
#    return fishmat
       
#if i == 0 : 
#                    inte1 = lambda k : float(basepk(k))**2*mp.cos(k*freqm)
#                    inte2 = lambda k : float(basepk(k))**2*mp.cos(k*freqp)
#                else : 
#                    inte1 = lambda k : basepk(k)*pkarr(k)*mp.cos(k*freqm)
#                    inte2 = lambda k : basepk(k)*pkarr(k)*mp.cos(k*freqp)
#                if j == l : 
#                    covm_w_der[i, j, l] = 1. - mp.quadosc(inte2, [0, mp.inf], zeros = lambda n : (n+0.5)*mp.pi/freqp)                
#                else : 
#                    covm_w_der[i, j, l] = mp.quadosc(inte1, [0, mp.inf], zeros = lambda n : (n+0.5)*mp.pi/freqm) \
#                                        - mp.quadosc(inte2, [0, mp.inf], zeros = lambda n : (n+0.5)*mp.pi/freqp)       