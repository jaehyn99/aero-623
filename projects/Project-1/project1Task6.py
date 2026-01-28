#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 17:50:06 2026

@author: curtis
"""

import numpy as np
import matplotlib.pyplot as plt
import project1Task4

#%%
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

def writegri(fname, V, E, B, Bname):
    Nn, dim = V.shape
    Ne = E.shape[0]

    with open(fname, "w") as f:
        # header
        f.write(f"{Nn} {Ne} {dim}\n")

        # vertices
        for v in V:
            f.write(" ".join(f"{x:.16e}" for x in v) + "\n")

        # boundaries
        NB = len(B)
        f.write(f"{NB}\n")
        for Bi, name in zip(B, Bname):
            Nb = Bi.shape[0]
            f.write(f"{Nb} 0 {name}\n")
            for row in Bi:
                f.write(" ".join(str(i + 1) for i in row) + "\n")

        # elements (single block)
        f.write(f"{Ne}\n")
        for row in E:
            f.write(" ".join(str(i + 1) for i in row) + "\n")

#-----------------------------------------------------------
def plotmesh(Mesh, fname):
    V = Mesh['vert']; E = Mesh['elem']; 
    f = plt.figure(figsize=(12,12))
    #plt.tripcolor(V[:,0], V[:,1], triangles=E)
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    #plt.plot(V[Mesh['bounds'][4][:, 0]][:, 0],V[Mesh['bounds'][4][:, 0]][:, 1], 'ro', markersize = 10)
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
    
    midpoint = {}
    B = Mesh['bounds']
    nF = np.size(B[0])//2
    refinedB = np.zeros([2*nF, 2])
    for jj in range(nF):
        if (B[0][jj][0] < B[0][jj][1]):
            key = (B[0][jj][0], B[0][jj][1])
            
        else:
            key = (B[0][jj][1], B[0][jj][0])
        
        x1 = V[B[0][jj][0], 0]; y1 = V[B[0][jj][0], 1]
        x2 = V[B[0][jj][1], 0]; y2 = V[B[0][jj][1], 1]
        
        xm = 0.5*(x1 + x2); ym = 0.5*(y1 + y2)
        midpoint[key] = np.size(V[:, 0])
        sS = 0
        sS, _ = project1Task4.gss(project1Task4.f2D, sS, project1Task4.sBreak, xm, ym, project1Task4.X, project1Task4.Y, project1Task4.S, project1Task4.dXi, project1Task4.dYi)
        xm = project1Task4.splineFun(sS, project1Task4.X, project1Task4.S, project1Task4.dXi); ym = project1Task4.splineFun(sS, project1Task4.Y, project1Task4.S, project1Task4.dYi)
        V = np.append(V, np.array([[xm, ym]]), axis = 0)

    nF = np.size(B[4])//2
    for jj in range(nF):
        if (B[4][jj][0] < B[4][jj][1]):
            key = (B[4][jj][0], B[4][jj][1])
            
        else:
            key = (B[4][jj][1], B[4][jj][0])
        
        x1 = V[B[4][jj][0], 0]; y1 = V[B[4][jj][0], 1]
        x2 = V[B[4][jj][1], 0]; y2 = V[B[4][jj][1], 1]
        
        xm = 0.5*(x1 + x2); ym = 0.5*(y1 + y2)
        midpoint[key] = np.size(V[:, 0])
        sS, _ = project1Task4.gss(project1Task4.f2D, sS, project1Task4.S[-1], xm, ym, project1Task4.X, project1Task4.Y + 18, project1Task4.S, project1Task4.dXi, project1Task4.dYi)
        xm = project1Task4.splineFun(sS, project1Task4.X, project1Task4.S, project1Task4.dXi); ym = project1Task4.splineFun(sS, project1Task4.Y + 18, project1Task4.S, project1Task4.dYi)
        V = np.append(V, np.array([[xm, ym]]), axis = 0)
        plt.plot(x1, y1, 'ko', markersize = 2)
        plt.plot(xm, ym, 'go', markersize = 3)
        plt.plot(x2, y2, 'ro')
    
    # Bisect all edges
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

M = readgri('/home/curtis/Documents/UMich/Courses/AEROSP623/project1/mesh_coarse.gri')
#Mesh = readgri('Z:/dotfiles/Documents/AEROSP623/project1/mesh_refined_2394.gri')
plotmesh(M, [])
rMesh = globalRefine(M)
plotmesh(rMesh, [])
#writegri('Z:/dotfiles/Documents/AEROSP623/project1/test2.gri', rMesh['vert'], rMesh['elem'], rMesh['bounds'], rMesh['Bname'])
#M2 = readgri('/home/curtis/Documents/UMich/Courses/AEROSP623/project1/test2.gri')
#plotmesh(M2, [])
#plotmesh(rMesh, []);
#rrMesh = globalRefine(rMesh)
#plotmesh(rrMesh, []);
#rrrMesh = globalRefine(rrMesh)
#plotmesh(rrrMesh, []);
