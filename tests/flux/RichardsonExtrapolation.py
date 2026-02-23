#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 22 16:02:02 2026

@author: curtis
"""

import os
import numpy as np
import matplotlib.pyplot as plt

#%%
base_dir = os.path.dirname(os.path.abspath(__file__))

cx = np.zeros([4])
cx[0] = 0.884120879
cx[1] = 0.844057964
cx[2] = 0.823268668
cx[3] = 0.814844963

cy = np.zeros([4])
cy[0] = 0.685169311
cy[1] = 0.739448345
cy[2] = 0.763917249
cy[3] = 0.767859289

cxP = np.log((cx[3] - cx[2])/(cx[2] - cx[1]))/np.log(2)
cyP = np.log((cy[3] - cy[2])/(cy[2] - cy[1]))/np.log(2)
cxRich = (2**cxP*cx[1] - cx[2])/(2**cxP - 1)
cyRich = (2**cyP*cy[1] - cy[2])/(2**cyP - 1)

dx = np.array([8, 4, 2, 1])
fig, axes = plt.subplots(2, 1, figsize=(6, 8), sharex=True)
axes[0].plot(dx, np.abs(cxRich - cx), 'ro', label = 'Simulation')
axes[0].plot(dx, dx**(-cxP)*np.abs(cxRich - cx)[-1], label = 'p = '+str(round(cxP, 3)))
axes[0].set_yscale('log'); axes[0].set_xscale('log')
axes[0].set_ylabel(r'|$c_{x\ error}$|')
axes[0].grid(); axes[0].legend()
axes[1].plot(dx, np.abs(cyRich - cy), 'ro')
axes[1].plot(dx, dx**(-cyP)*np.abs(cyRich - cy)[-1], label = 'p = '+str(round(cyP, 3)))
axes[1].set_yscale('log'); axes[0].set_xscale('log')
axes[1].set_ylabel(r'|$c_{y\ error}$|')
axes[1].set_xlabel(r'$\Delta x / \Delta x_{finest}$')
axes[1].grid(); axes[1].legend(); 
filePath = os.path.join(base_dir, "pressureCoeffConvergence.png")
plt.savefig(filePath, dpi=250)