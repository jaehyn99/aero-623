#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 16 08:37:55 2026

@author: curtis
"""

import os
import numpy as np
from scipy.optimize import bisect
import matplotlib.pyplot as plt

#%% Compute the analytical solution to the Sod shock tube problem

def res(p, pL, rhoL, pR, rhoR, gamma, aL, AR, BR):
    return 2*aL*((p/pL)**(gm1/(2*gamma)) - 1)/gm1 + (p - pR)*np.sqrt(AR/(p + BR))

pL = 1; rhoL = 1; pR = 0.1; rhoR = 0.125; gamma = 1.4
gm1 = gamma - 1; gp1 = gamma + 1
aL = np.sqrt(gamma*pL/rhoL); gm1 = gamma - 1
AR = 2/(gp1*rhoR); BR = gm1*pR/gp1
aL = np.sqrt(gamma*pL/rhoL)
pS = bisect(res, 0.1, 1, args = (pL, rhoL, pR, rhoR, gamma, aL, AR, BR))
uS = -0.5*(2*aL*((pS/pL)**(gm1/(2*gamma)) - 1)/gm1 - (pS - pR)*np.sqrt(AR/(pS + BR)))
aSL = aL*(pS/pL)**(gm1/(2*gamma))
SHL = -aL; STL = uS - aSL
QR = np.sqrt((pS + BR)/AR)
SR = QR/rhoR
x = np.linspace(-0.5, 0.5, 1000)
t = 0.25
rho = np.zeros_like(x)
u = np.zeros_like(x)
p = np.zeros_like(x)
for ii in range(np.size(x)):
    if (x[ii]/t < SHL):
        rho[ii], u[ii], p[ii] = rhoL, 0, pL
    
    elif (x[ii]/t < STL):
        rhoLFan = rhoL*(2/gp1 + gm1/(gp1*aL)*(-x[ii]/t))**(2/gm1)
        uLFan = 2/gp1*(aL + x[ii]/t)
        pLFan = pL*(2/gp1 + gm1/(gp1*aL)*(-x[ii]/t))**(2*gamma/gm1)
        rho[ii], u[ii], p[ii] = rhoLFan, uLFan, pLFan
        
    elif (x[ii]/t < uS):
        rhoSL = rhoL*(pS/pL)**(1/gamma)
        rho[ii], u[ii], p[ii] = rhoSL, uS, pS
        
    elif (x[ii]/t < SR):
        rhoSR = rhoR*((pS/pR + gm1/gp1)/(gm1/gp1*pS/pR + 1))
        rho[ii], u[ii], p[ii] = rhoSR, uS, pS
        
    else:
        rho[ii], u[ii], p[ii] = rhoR, 0, pR
        
#%% Plot solutions from approximate Riemann solvers over the theoretical solution

base_dir = os.path.dirname(os.path.abspath(__file__))

filePath = os.path.join(base_dir, "sodHLLE.csv")
xHLLE = np.loadtxt(filePath, usecols=(0,), delimiter = ',')
rhoHLLE = np.loadtxt(filePath, usecols=(1,), delimiter = ',')
uHLLE = np.loadtxt(filePath, usecols=(2,), delimiter = ',')/rhoHLLE
pHLLE = (np.loadtxt(filePath, usecols=(4,), delimiter = ',') - 0.5*rhoHLLE*uHLLE**2)*gm1

filePath = os.path.join(base_dir, "sodRoe.csv")
xRoe = np.loadtxt(filePath, usecols=(0,), delimiter = ',')
rhoRoe = np.loadtxt(filePath, usecols=(1,), delimiter = ',')
uRoe = np.loadtxt(filePath, usecols=(2,), delimiter = ',')/rhoRoe
pRoe = (np.loadtxt(filePath, usecols=(4,), delimiter = ',') - 0.5*rhoRoe*uRoe**2)*gm1

fig, axes = plt.subplots(3, 1, figsize=(6, 8), sharex=True)
axes[0].plot(xHLLE, rhoHLLE, 'r-', label = 'HLLE')
axes[0].plot(xRoe, rhoRoe, 'b-', label = 'Roe')
axes[0].plot(x, rho, 'k-', label = 'Theoretical')
axes[0].set_ylabel(r"$\rho$")
axes[0].grid()

axes[1].plot(xHLLE, uHLLE, 'r-', label = 'HLLE')
axes[1].plot(xRoe, uRoe, 'b-', label = 'Roe')
axes[1].plot(x, u, 'k-', label = 'Theoretical')
axes[1].set_ylabel("u")
axes[1].legend()
axes[1].grid()

axes[2].plot(xHLLE, pHLLE, 'r-', label = 'HLLE')
axes[2].plot(xRoe, pRoe, 'b-', label = 'Roe')
axes[2].plot(x, p, 'k-', label = 'Theoretical')
axes[2].set_ylabel("p")
axes[2].set_xlabel("x")
axes[2].grid()

plt.tight_layout()
filePath = os.path.join(base_dir, "sodComparison.png")
plt.savefig(filePath, dpi=250)
