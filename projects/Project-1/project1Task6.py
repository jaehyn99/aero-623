#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 17:50:06 2026

@author: curtis
"""

import numpy as np
import matplotlib.pyplot as plt

#-----------------------------------------------------------
def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    V = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    B = []; Bname = []
    for i in range(NB):
        s = f.readline().split(); Nb = int(s[0]); Bname.append(s[2])
        Bi = np.array([[int(s)-1 for s in f.readline().split()] for n in range(Nb)])
        B.append(Bi)
    # read elements
    Ne0 = 0; E = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0==0) else np.concatenate((E,Ei), axis=0)
        Ne0 += ne
    f.close()
    Mesh = {'vert':V, 'elem':E, 'bounds':B, 'Bname':Bname }
    return Mesh

#-----------------------------------------------------------
def plotmesh(Mesh, fname):
    V = Mesh['vert']; E = Mesh['elem']; 
    f = plt.figure(figsize=(12,12))
    #plt.tripcolor(V[:,0], V[:,1], triangles=E)
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    dosave = not not fname
    plt.axis('equal')
    plt.tick_params(axis='both', labelsize=12)
    f.tight_layout(); plt.show(block=(not dosave))
    if (dosave): plt.savefig(fname)
    plt.close(f)

#-----------------------------------------------------------

def globalRefine(Mesh):
    V = Mesh['vert']
    E = Mesh['elem']
    nE = np.size(E[:, 0])
    nB= len(Mesh['bounds'])
    # Bisect all edges
    midpoint = {}
    for ii in range(nE):
        for jj in range(3):
            if (E[ii, jj - 1] < E[ii, jj]):
                key = (E[ii, jj - 1], E[ii, jj])
                
            else:
                key = (E[ii, jj], E[ii, jj - 1])
                
            if key in midpoint:
                continue
            
            x1 = V[E[ii, jj - 1], 0]; y1 = V[E[ii, jj - 1], 1]
            x2 = V[E[ii, jj], 0]; y2 = V[E[ii, jj], 1]
            
            xm = 0.5*(x1 + x2); ym = 0.5*(y1 + y2)
            midpoint[key] = np.size(V[:, 0])
            V = np.append(V, np.array([[xm, ym]]), axis = 0)
    
    # Define refined mesh with new verticies
    refinedE = np.zeros([4*nE, 3])
    key = np.zeros([3, 2])
    for ii in range(nE):
        for jj in range(3):
            if (E[ii, jj - 2] < E[ii, jj - 1]):
                key[jj - 1, :] = (E[ii, jj - 2], E[ii, jj - 1])
                
            else:
                key[jj - 1, :] = (E[ii, jj - 1], E[ii, jj - 2])
        
        for jj in range(3):    
            refinedE[4*ii + jj, -1] = midpoint[tuple(key[jj - 1, :])]
            refinedE[4*ii + jj, 0] = E[ii, jj - 1]
            refinedE[4*ii + jj, 1] = midpoint[tuple(key[jj, :])]
            refinedE[4*ii + 3, jj] = midpoint[tuple(key[jj, :])]
    
    # Update bounds
    B = Mesh['bounds']
    for ii in range(nB):
        nF = np.size(B[ii])//2
        refinedB = np.zeros([2*nF, 2])
        for jj in range(nF):
            if (B[ii][jj][0] < B[ii][jj][1]):
                key = (B[ii][jj][0], B[ii][jj][1])
                
            else:
                key = (B[ii][jj][1], B[ii][jj][0])
                
            refinedB[2*jj] = np.array([B[ii][jj][0], midpoint[key]])
            refinedB[2*jj + 1] = np.array([midpoint[key], B[ii][jj][1]])
        
        B[ii] = refinedB.astype('int')
    
    refinedMesh = {'vert':V, 'elem':refinedE.astype('int'), 'bounds':B, 'Bname':Mesh['Bname']}
    return refinedMesh

Mesh = readgri('/home/curtis/Documents/UMich/Courses/AEROSP623/project1/test.gri')
refinedMesh = globalRefine(Mesh)
plotmesh(refinedMesh, []);
