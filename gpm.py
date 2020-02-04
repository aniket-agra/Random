#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:03:56 2019

@author: aagrawal
"""

import george
from george import kernels
import emcee
from george.modeling import Model
import corner

correls_d = []
for i in range(1,25):
    for j in range(1,4):
        d=np.loadtxt('%02d/2pcf/all/1gpc/xi_cra0_grd0_a__los%d_z014' %(i,j)).T
        correls_d.append(d[2,]/d[1,]-1)
correls_d = np.array(correls_d)

correls_p = []
for i in range(1,25):
    for j in range(1,4):
        d=np.loadtxt('%02d/2pcf/all/1gpc/xi_cra1_grd1_a__los%d_z014' %(i,j)).T
        correls_p.append(d[2,]/d[1,])
correls_p = np.array(correls_p)

correls_v = []
for i in range(1,25):
    for j in range(1,4):
        d=np.loadtxt('%02d/2pcf/all/1gpc/xi_cra1_grd1_v__los%d_z014' %(i,j)).T
        correls_v.append(d[2,]/d[1,])
correls_v = np.array(correls_v)

x = d[0,]

x = x[:40]
correls_d = correls_d[:,:40]
correls_p = correls_p[:,:40]
correls_v = correls_v[:,:40]
print('separations:', x, 'in Mpc/h')

data = np.concatenate((correls_d,correls_p,correls_v),axis=1)

cov = np.cov(data.T)
R_coeff = np.corrcoef(data.T)
plt.figure(figsize=(10,10))
plt.title('correlation coefficient $C_{ij}\,/\,\sqrt{C_{ii}C_{jj}}$')
plt.imshow(R_coeff)

ystd_d = np.std(np.log10(correls_d),axis=0,ddof=1)
ystd_p = np.std(np.log10(correls_p),axis=0,ddof=1)
ystd_v = np.std(np.log10(correls_v),axis=0,ddof=1)

class LinModel(Model):
    parameter_names = ("zeroth", "first")

    def get_value(self, t):
        return self.first * t.flatten()+self.zeroth
    
def lnprob(p):
    gp.set_parameter_vector(p)
    return gp.log_likelihood(y, quiet=True)

i, j = 1, 1 # consider the 1st realization, first projection direction

y = np.log10(correls_d[3*(i-1)+j-1,]) # output
t = np.log10(x)                # input

inipara = dict(zeroth=0.,first=-1.) # initial guess of the linear dependence
kernel = np.var(y)*kernels.ExpSquaredKernel(0.01) # gaussian kernel assumed
gp = george.GP(kernel,mean=LinModel(**inipara))

gp.compute(t,ystd_d) # give it the uncertainty of the data. 
                                        # Off-diagonal component is not supported currently, 
                                        # but GP kernel is expected to capture the covariance between data points.

initial = gp.get_parameter_vector() # the full parameter vector; kernel parameters + linear model parameters

# perform MCMC
ndim, nwalkers = len(initial), 32
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob)

print("Running first burn-in...")
p0 = initial + 1e-6 * np.random.randn(nwalkers, ndim)
p0, lp, _ = sampler.run_mcmc(p0, 2000)

print("Running second burn-in...")
p0 = p0[np.argmax(lp)] + 1e-6 * np.random.randn(nwalkers, ndim)
sampler.reset()
p0, _, _ = sampler.run_mcmc(p0, 1000)
sampler.reset()

print("Running production...")
sampler.run_mcmc(p0, 1000);

short_param_names = ('y-intercept','slope','$\log(A)$','$\log(\sigma^2)$')
plt.figure()
corner.corner(sampler.flatchain,labels=short_param_names)
plt.show()

pbest = sampler.flatchain[np.argmax(sampler.flatlnprobability),:]
print 'best-fit params:', pbest
gp.set_parameter_vector(pbest)

ts = np.linspace(-1.,1.,500)
ys,cov = gp.predict(y,ts,return_cov=True)
err = np.sqrt(diag(cov))

plt.figure(figsize=(10,6))
plt.plot(ts,ys)
plt.plot(ts,gp.mean.get_value(ts), label='unconstrained mean function (assumed to be linear)')
plt.fill_between(ts,(ys+err),(ys-err),color='g',alpha=0.5, label='constrained mean with $1-\sigma$ uncertainty')
plt.errorbar(t,y,yerr=ystd_d,fmt='ok',capsize=4, label='data to be learned')
plt.xlabel('$\log_{10}(x\,/\,[h^{-1}\mathrm{Mpc}])$')
plt.ylabel('$\log_{10}(\\xi_{gc})$')
plt.legend(loc='upper right',fontsize=14)
plt.show()

logslope = np.gradient(ys)/np.gradient(ts)
R_sp = 10**(ts[np.argmin(logslope)])
minslope = logslope[np.argmin(logslope)]
plt.figure(figsize=(10,6))
plt.plot(ts,logslope)
plt.xlabel('$\log_{10}(x\,/\,[h^{-1}\mathrm{Mpc}])$')
plt.ylabel('$\mathrm{d}\log_{10}(\\xi_{gc})\, / \, \mathrm{d} \log_{10}(x\,/\,[h^{-1}\mathrm{Mpc}])$')
print('minimum slope:', minslope, 'at', R_sp, '[Mpc/h]')

