#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 19:46:25 2019

@author: aagrawal
"""
import GPy
import numpy as np
from matplotlib import pyplot as plt

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['xtick.major.size'] = 5
plt.rcParams['ytick.major.size'] = 5
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

correls_d = []
d = np.loadtxt('xi_cra0_mrd0_a_0512_los1_z014')
#correls_d = np.array(correls_d)

x = d[:,0]

x = x[:40]

print('separations:', x, 'in Mpc/h')

t = np.log10(x)
y = np.log10(d[0:40,1]).T
#plt.plot(t,y[0,],'o')
#plt.plot(t,y[1,],'^')
#plt.plot(t,y[2,],'v')

K1 = GPy.kern.RBF(1)
K2 = GPy.kern.Bias(1) + GPy.kern.Linear(1)
B1 = GPy.kern.Coregionalize(1, output_dim=1)
B2 = GPy.kern.Coregionalize(1, output_dim=1)
kernel = K1**B1 + K2**B2

model = GPy.models.GPCoregionalizedRegression(X_list=[t[:,None]], Y_list=[y[:,None]], kernel=kernel)
model.optimize()
model.plot(plot_limits=[-1, 1], fixed_inputs=[(1, 0)],which_data_rows=slice(0,40))
#model.plot(plot_limits=[-1, 1], fixed_inputs=[(1, 1)],which_data_rows=slice(40,80))
#model.plot(plot_limits=[-1, 1], fixed_inputs=[(1, 2)],which_data_rows=slice(80,120))

ts = np.linspace(-1,1,200)
newX = ts[:,None]
newX = np.hstack([newX,np.ones_like(newX)])
noise_dict = {'output_index':newX[:,1:].astype(int)}

mu, var = model.predict(newX,Y_metadata=noise_dict)
plt.figure(figsize=(10,6))
plt.errorbar(ts,mu,yerr=var,fmt='o',capsize=4)
