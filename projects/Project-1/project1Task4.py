#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 12 18:16:34 2026

@author: curtis
"""

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#%% 2D Spline

def tridiag(A, B, C, D):
    n = np.size(A)

    CP = np.zeros(n-1)
    DP = np.zeros(n)
    x  = np.zeros(n)

    CP[0] = C[0] / A[0]
    DP[0] = D[0] / A[0]

    for i in range(1, n-1):
        denom = A[i] - B[i-1] * CP[i-1]
        CP[i] = C[i] / denom
        DP[i] = (D[i] - B[i-1] * DP[i-1]) / denom

    denom = A[n-1] - B[n-2] * CP[n-2]
    DP[n-1] = (D[n-1] - B[n-2] * DP[n-2]) / denom

    x[n-1] = DP[n-1]
    for i in range(n-2, -1, -1):
        x[i] = DP[i] - CP[i] * x[i+1]

    return x

def splineFit(X, S):
    A = np.zeros(np.size(X))
    B = np.zeros(np.size(X) - 1)
    C = np.zeros(np.size(X) - 1)
    D = np.zeros(np.size(X))

    DS = np.zeros(np.size(S) - 1)
    for ii in range(np.size(DS)):
        DS[ii] = S[ii + 1] - S[ii]
        
    A[0] = 2*DS[0]
    C[0] = DS[0]
    D[0] = 3*(X[1] - X[0])
    for ii in range(1, np.size(A) - 1):
        B[ii - 1] = DS[ii]
        A[ii] = 2*(DS[ii - 1] + DS[ii])
        C[ii] = DS[ii - 1]
        D[ii] = 3*((X[ii] - X[ii - 1])/DS[ii - 1]*DS[ii] + (X[ii + 1] - X[ii])/DS[ii]*DS[ii - 1])
        
    B[-1] = DS[-2]
    A[-1] = 2*DS[-2]
    D[-1] = 3*(X[-1] - X[-2])
    dX = tridiag(A, B, C, D)
    
    return dX 

def splineFun(s, Xi, Si, dXi):
    lo, hi = 0, len(Si)
    while lo < hi:
        ind = (lo + hi) // 2
        if Si[ind] < s:
            lo = ind + 1
        else:
            hi = ind
    
    ind = (lo + hi) // 2 - 1        
    DSi = Si[ind + 1] - Si[ind]
    t = (s - Si[ind])/DSi
    xi0P = (dXi[ind] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    xi1P = (dXi[ind + 1] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    x = (1 - t)*Xi[ind] + t*Xi[ind + 1] + (t - t**2)*((1 - t)*xi0P - t*xi1P)
   
    return x

def diffSplineFun(s, Xi, Si, dXi):
    lo, hi = 0, len(Si)
    while lo < hi:
        ind = (lo + hi) // 2
        if Si[ind] < s:
            lo = ind + 1
        else:
            hi = ind
    
    ind = (lo + hi) // 2 - 1
    DSi = Si[ind + 1] - Si[ind]
    t = (s - Si[ind])/DSi
    xi0P = (dXi[ind] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    xi1P = (dXi[ind + 1] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    dx = (Xi[ind + 1] - Xi[ind] + (1 - 2*t)*((1 - t)*xi0P - t*xi1P) - (t - t**2)*(xi0P + xi1P))
    return dx/DSi

def diff2SplineFun(s, Xi, Si, dXi):
    lo, hi = 0, len(Si)
    while lo < hi:
        ind = (lo + hi) // 2
        if Si[ind] < s:
            lo = ind + 1
        else:
            hi = ind
    
    ind = (lo + hi) // 2 - 1
    DSi = Si[ind + 1] - Si[ind]
    t = (s - Si[ind])/DSi
    xi0P = (dXi[ind] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    xi1P = (dXi[ind + 1] - (Xi[ind + 1] - Xi[ind])/DSi)*DSi
    d2x = (6*t-4)*xi0P + (6*t-2)*xi1P
    return d2x/(DSi*DSi)

def sSimps(Xi, Xip1, Yi, Yip1, Si, Sip1, dXi, dXip1, dYi, dYip1):
    xi0P = (dXi - (Xip1 - Xi)/(Sip1 - Si))*(Sip1 - Si)
    xi1P = (dXip1 - (Xip1 - Xi)/(Sip1 - Si))*(Sip1 - Si)
    yi0P = (dYi - (Yip1 - Yi)/(Sip1 - Si))*(Sip1 - Si)
    yi1P = (dYip1 - (Yip1 - Yi)/(Sip1 - Si))*(Sip1 - Si)
    dxi0 = Xip1 - Xi + xi0P
    dyi0 = Yip1 - Yi + yi0P
    f0 = np.sqrt(dxi0**2 + dyi0**2)
    dxi1 = Xip1 - Xi - xi0P/4 - xi1P/4
    dyi1 = Yip1 - Yi - yi0P/4 - yi1P/4
    f1 = np.sqrt(dxi1**2 + dyi1**2)
    dxi2 = Xip1 - Xi + xi1P
    dyi2 = Yip1 - Yi + yi1P
    f2 = np.sqrt(dxi2**2 + dyi2**2)
    return (f0 + 4*f1 + f2)/6

def splineFit2D(X, Y, eps = 1E-12):
    sizeS = np.size(X)
    S = np.zeros(sizeS)
    for ii in range(1, sizeS):
        S[ii] = S[ii - 1] + np.sqrt((X[ii] - X[ii - 1])**2 + (Y[ii] - Y[ii - 1])**2)
        
    splineErr = 1
    while (splineErr > eps):
        dXi = splineFit(X, S)
        dYi = splineFit(Y, S)
        splineErr = 0
        trueS = np.zeros(sizeS)
        for ii in range(1, sizeS):
            trueS[ii] = trueS[ii - 1] + sSimps(X[ii - 1], X[ii], Y[ii - 1], Y[ii], S[ii - 1], S[ii], dXi[ii - 1], dXi[ii], dYi[ii - 1], dYi[ii])
            splineErr += np.abs(trueS[ii] - S[ii])
            
        S = trueS
        splineErr /= S[-1]
        
    return S, dXi, dYi

def sBreakRes(s, xU, X, S, dXi):
    return xU[-1] - splineFun(s, X, S, dXi)

gr = (np.sqrt(5.0) - 1.0) / 2.0

def f2D(s, p, q, Xi, Yi, Si, dXi, dYi):
    return (splineFun(s, Xi, Si, dXi) - p)**2 + (splineFun(s, Yi, Si, dYi) - q)**2

def gss(f, a, b, p, q, Xi, Yi, Si, dXi, dYi, tol=1e-14, max_iter=200):
    c = b - gr * (b - a)
    d = a + gr * (b - a)
    fc = f(c, p, q, Xi, Yi, Si, dXi, dYi)
    fd = f(d, p, q, Xi, Yi, Si, dXi, dYi)

    for _ in range(max_iter):
        if abs(b - a) < tol:
            break

        if fc < fd:
            b = d
            d = c
            fd = fc
            c = b - gr * (b - a)
            fc = f(c, p, q, Xi, Yi, Si, dXi, dYi)
        else:
            a = c
            c = d
            fc = fd
            d = a + gr * (b - a)
            fd = f(d, p, q, Xi, Yi, Si, dXi, dYi)

    x_best = 0.5 * (a + b)
    return x_best, f(x_best, p, q, Xi, Yi, Si, dXi, dYi)

xU = np.loadtxt('bladeupper.txt', usecols=(0,))
yU = np.loadtxt('bladeupper.txt', usecols=(1,))
# xL = np.loadtxt('/home/curtis/Documents/UMich/Courses/AEROSP623/project1/bladelower.txt', usecols=(0,))
# yL = np.loadtxt('/home/curtis/Documents/UMich/Courses/AEROSP623/project1/bladelower.txt', usecols=(1,))

X = xU # np.append(xU, xL[1:])
Y = yU # np.append(yU, yL[1:])

S, dXi, dYi = splineFit2D(X, Y)

xPlot = np.zeros(200)
yPlot = np.zeros(200)
sPlot = np.linspace(S[0], S[-1], 200)

for ii in range(200):
   xPlot[ii] = splineFun(sPlot[ii], X, S, dXi)
   yPlot[ii] = splineFun(sPlot[ii], Y, S, dYi)
plt.plot(xPlot, yPlot)

# range of random points

print(S[0], S[-1])

px = [3.1677, -4.7468, -4.2695, -4.7723, -3.8026, -8.42, 1.0205, 3.9801, 5.7149, -2.2342]
py = [5.3846, 3.9693, -4.7231, 3.0665, 5.5248, -1.9938, -5.4757, -7.7584, -11.3292, -11.6956]

for i in range(10):
    p = [px[i], py[i]]
    s, _ = gss(f2D, 0, S[-1], p[0], p[1], X, Y, S, dXi, dYi)

    r = np.array([splineFun(s, X, S, dXi), splineFun(s, Y, S, dYi)])
    plt.scatter(p[0], p[1], color="black")
    plt.plot((p[0], r[0]), (p[1], r[1]), color="red")

    dr = np.array([diffSplineFun(s, X, S, dXi), diffSplineFun(s, Y, S, dYi)])
    d2r = np.array([diff2SplineFun(s, X, S, dXi), diff2SplineFun(s, Y, S, dYi)])

    dLds = 2*np.dot(r-p, dr)
    dL2ds = 2*np.dot(dr, dr) + 2*np.dot(r-p, d2r)

    print(p, s, np.linalg.norm(r-p), dLds, dL2ds)

plt.xlabel("x")
plt.ylabel("y")
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig("spline_projection_test.png")

# sBreak = fsolve(sBreakRes, 30, args = (xU, X, S, dXi))[0]
# sPlot = np.linspace(S[0], sBreak, 100)
# for ii in range(100):
#     xPlot[ii] = splineFun(sPlot[ii], X, S, dXi)
#     yPlot[ii] = splineFun(sPlot[ii], Y, S, dYi)

# plt.plot(xPlot, yPlot, 'ro', markersize = 2)
# sPlot = np.linspace(sBreak, S[-1], 100)
# for ii in range(100):
#     xPlot[ii] = splineFun(sPlot[ii], X, S, dXi)
#     yPlot[ii] = splineFun(sPlot[ii], Y + 18, S, dYi)

# plt.plot(xPlot, yPlot, 'ro', markersize = 2)
# plt.plot(xU, yU, 'ko', markersize = 1)
# plt.plot(xL, yL + 18, 'ko', markersize = 1)
# ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')
# plt.show()
# #%% 2D Spline Projection
# gr = (np.sqrt(5.0) - 1.0) / 2.0

# p = 15; q = -5

# plt.plot(p, q, 'ro'); 
# sPlot = np.linspace(S[0], sBreak, 100)
# for ii in range(100):
#     xPlot[ii] = splineFun(sPlot[ii], X, S, dXi)
#     yPlot[ii] = splineFun(sPlot[ii], Y, S, dYi)

# plt.plot(xPlot, yPlot, 'ro', markersize = 2)
# sPlot = np.linspace(sBreak, S[-1], 100)
# for ii in range(100):
#     xPlot[ii] = splineFun(sPlot[ii], X, S, dXi)
#     yPlot[ii] = splineFun(sPlot[ii], Y + 18, S, dYi)

# plt.plot(xPlot, yPlot, 'ro', markersize = 2)

# s1, _ = gss(f2D, 0, sBreak, p, q, X, Y, S, dXi, dYi)
# s2, _ = gss(f2D, sBreak, S[-1], p, q, X, Y + 18, S, dXi, dYi)

# if f2D(s1, p, q, X, Y, S, dXi, dYi) <= f2D(s2, p, q, X, Y + 18, S, dXi, dYi):
#     sS = s1
#     plt.plot([p, splineFun(sS, X, S, dXi)], [q, splineFun(sS, Y, S, dYi)])

# else:
#     sS = s2
#     plt.plot([p, splineFun(sS, X, S, dXi)], [q, splineFun(sS, Y + 18, S, dYi)])

# ax = plt.gca()
# ax.set_aspect('equal', adjustable='box')
# plt.show()
