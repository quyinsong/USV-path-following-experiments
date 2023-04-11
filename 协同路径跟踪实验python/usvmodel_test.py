# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 16:31:33 2022

@author: qys
"""

import numpy as np
import matplotlib.pyplot as plt
from usvmodel import USV1

ts = 0.02
tfinal = 40
Ns = tfinal/ts
Tc = np.array([[100],[50]])
xout = np.zeros((int(Ns),13))

for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        nu = np.array([[0],[0],[0]])
        eta = np.array([[0],[0],[0]])
        Tf = np.array([[0],[0]])
    
    # disturbancs and currents
    tw = np.array([[0],[0],[0]])
    vc = np.array([[0],[0],[0]])
    
    # USV
    usv1 = USV1()
    nu,eta,Tf,fu,fr = usv1.state_update(nu,eta,Tf,Tc,vc,tw,ts)
    
    # store time series
    xout[i,:] = np.hstack((np.array([t]),nu[0],nu[1],nu[2],eta[0],eta[1],eta[2],Tc[0],Tc[1],
                           Tf[0],Tf[1],fu[0],fr[0]))
    
t = xout[:,0]
u = xout[:,1]
v = xout[:,2]
r = xout[:,3]
x = xout[:,4]
y = xout[:,5]
psi = xout[:,6]
Tc1 = xout[:,7]
Tc2 = xout[:,8]
Tf1 = xout[:,9]
Tf2 = xout[:,10]
fu = xout[:,11]
fr = xout[:,12]

# plots
font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12,
}

fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(y,x)
plt.xlabel("east (m)",font)
plt.ylabel("north (m)",font)
plt.show()

fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(t,Tc1)
plt.plot(t,Tc2)
plt.xlabel("time (s)")
plt.ylabel("thrust (N)")
plt.show()

fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,u)
plt.ylabel("u (m/s)",font,labelpad=18)

plt.subplot(312)
axarr[1].plot(t,v)
plt.ylabel("v (m/s)",font,labelpad=6)

plt.subplot(313)
axarr[2].plot(t,r)
plt.ylabel("r (rad/s)",font,labelpad=18)

plt.xlabel("time (s)", font)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.2)





















        
    
    