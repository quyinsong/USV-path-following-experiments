# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 09:59:02 2022

@author: qys
"""

import numpy as np
import matplotlib.pyplot as plt
from usvmodel import USV1
from nncontroller import NNcontroller1

# initialize
ts = 0.02
tfinal = 40
Ns = tfinal/ts

Nxout = 2
Nusv1 = 14

xout = np.zeros((int(Ns),Nxout))
usv1_out = np.zeros((int(Ns),Nusv1))

# path information
def fxy1(x,y):

    f = -2*x+y+32
    yd = 2*x-32
    fx = -2
    fy = 1
    ffxx = 0
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy,yd
    
def fxy2(x,y):
    f = -1*x+y
    yd = x
    fx = -1
    fy = 1
    ffxx = 0
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy,yd
    
def fxy3(x,y):
    f = -x+2*y
    yd = x/2
    fx = -1
    fy = 2
    ffxx = 0
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy,yd

def fxy4(x,y):
    Ro = 10.
    f = x**2+y**2-Ro**2
    yd = 0
    fx = 2*x
    fy = 2*y
    ffxx = 2
    ffyy = 2
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy,yd
    
path_info = {'p1':fxy1,
             'p2':fxy2,
             'p3':fxy3,
             'p4':fxy4}
    
def path(x,y,path_num):
    return path_info.get(path_num)(x,y)

# simulation start
for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        usv1_nu = np.array([[0],[0],[0]]) # usv1-usv3 states initialization
        usv1_eta = np.array([[0],[-20],[0]])
        usv1_Tf = np.array([[0],[0]])
        usv1_Tc = np.array([[50],[50]])
        usv2_nu = np.array([[0],[0],[0]])
        uc = np.array([0.8])
    
    # disturbancs and currents
    tw = np.array([[0],[0],[0]])
    vc = np.array([[0],[0],[0]])
    
    # instantiate USV object
    usv1 = USV1()
    usv1_nu,usv1_eta,usv1_nu_dot, usv1_eta_dot, usv1_Tf,usv1_fu,\
    usv1_fr = usv1.state_update1(usv1_nu,usv1_eta,usv1_Tf,usv1_Tc,vc,tw,ts)
    
    # single guidance law
    if t == 0:
        k1 = 0.2
        k2 = 0.8
        
    f,fx,fy,ffxx,ffyy,ffxy,yd = path(usv1_eta[0],usv1_eta[1],'p1')
    fxy_norm = np.linalg.norm(np.array([fx,fy]))
    psi = usv1_eta[2]
    thetad = -np.arctan(fx/fy)
    v = usv1_nu[1]
    u = usv1_nu[0]
    betad = np.arctan(v/uc)
    thetae = psi - thetad+betad
    f_dot = fxy_norm*usv1_nu[0]*np.sin(thetae)
    thetad_dot = ((fx*ffxy-fy*ffxx)*u*np.cos(psi)+\
                  (fx*ffyy-fy*ffxy)*u*np.sin(psi))/(fxy_norm**2)
    Ud = np.sqrt(uc**2+v**2)
    betad_dot = uc*usv1_nu_dot[1]/(Ud**2)
    s = thetae+np.arctan(k1*f)
    rc = thetad_dot-betad_dot-k1*f_dot/(1+(k1*f)**2)-k2*np.sign(s)*(np.abs(s)**0.9)
    
    # NN controller
    usv1_wd = np.array([uc,rc])

    usv1_Delta_T = usv1_Tc-usv1_Tf

    if t==0:
        B_usv = np.array([[1/50.5,1/50.5],[0.26/17.21,-0.26/17.21]])
        usv1_wdf = usv1_wd
        usv1_Lambda = np.dot(B_usv,usv1_Delta_T)
        usv1_What = np.zeros((9,2))
        
    nnc1 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,2]),A_aux=np.diag([1,0.5]))
        
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)

    usv1_X = np.array([usv1_nu[0],usv1_nu[1],usv1_nu[0]*usv1_nu[2]])

    usv1_w = np.array([usv1_nu[0],usv1_nu[2]])
     
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)
    
    usv1_ew = usv1_w-usv1_wd-usv1_Lambda
    
    usv1_wdf, usv1_wdf_dot = nnc1.TD(usv1_wd,usv1_wdf,ts)
    
    usv1_Fhat, usv1_What = nnc1.estimator(usv1_X,usv1_What,usv1_ew,ts)
    
    usv1_Tc = nnc1.control_law(usv1_ew,usv1_Fhat,usv1_wdf_dot,usv1_Lambda)
       
    # store time series
    xout[i,:] = np.hstack((np.array([t]),yd))
    usv1_out[i,:] = np.hstack((usv1_nu[0],usv1_nu[1],usv1_nu[2],\
                               usv1_eta[0],usv1_eta[1],usv1_eta[2],\
                               usv1_Tc[0],usv1_Tc[1],usv1_Tf[0],\
                               usv1_Tf[1],usv1_fu[0],usv1_fr[0],\
                               usv1_Fhat[0],usv1_Fhat[1]))
    
    
    
t = xout[:,0]
yd = xout[:,1]

usv1_u = usv1_out[:,0]
usv1_v = usv1_out[:,1]
usv1_r = usv1_out[:,2]
usv1_x = usv1_out[:,3]
usv1_y = usv1_out[:,4]
usv1_psi = usv1_out[:,5]
usv1_Tc1 = usv1_out[:,6]
usv1_Tc2 = usv1_out[:,7]
usv1_Tf1 = usv1_out[:,8]
usv1_Tf2 = usv1_out[:,9]
usv1_fu = usv1_out[:,10]
usv1_fr = usv1_out[:,11]
usv1_fuhat = usv1_out[:,12]
usv1_frhat = usv1_out[:,13]


# plots
font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12,
}

# north-east map
fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(usv1_y,usv1_x,'r-')
plt.plot(yd,usv1_x,'b--')
plt.xlabel("east (m)",font)
plt.ylabel("north (m)",font)

plt.savefig("NNtest_map.svg",dpi = 600, format = "svg",bbox_inches = "tight")

# thrust
fig = plt.figure(figsize=(6,3),dpi = 600)

plt.plot(t,usv1_Tc1,'r-',label = "Tc1")
plt.plot(t,usv1_Tc2,'b-',label = "Tc2")
plt.plot(t,usv1_Tf1,'g--',label = "T1")
plt.plot(t,usv1_Tf2,'m--',label = "T2")

plt.xlabel("time (s)",font,labelpad=18)
plt.ylabel("Thrust (N)",font,labelpad=18)
plt.legend(bbox_to_anchor = (0.8,1)) # legend position

plt.savefig("NNtest_thrust.svg",dpi = 600, format = "svg",bbox_inches = "tight")

# disturbances estimate
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(2,1,sharex = True)

plt.subplot(211)
axarr[0].plot(t,usv1_fu,'r-',label = "$f_{u}$")
axarr[0].plot(t,usv1_fuhat,'b--',label="$\hat{f}_{u}$")
plt.ylabel("fu (N)",font,labelpad=18)
plt.legend(loc="upper right")

plt.subplot(212)
axarr[1].plot(t,usv1_fr,'r-',label = "$f_{r}$")
axarr[1].plot(t,usv1_frhat,'b--',label = "$\hat{f}_{r}$")
plt.legend(loc="upper right")

plt.xlabel("time (s)",font,labelpad=18)
plt.ylabel("fr (Nm)",font,labelpad=18)

plt.savefig("NNtest_estimation.svg",dpi=600,format = "svg",bbox_inches = "tight")

# velocities
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,usv1_u)
plt.ylabel("u (m/s)",font,labelpad=18)

plt.subplot(312)
axarr[1].plot(t,usv1_v)
plt.ylabel("v (m/s)",font,labelpad=6)

plt.subplot(313)
axarr[2].plot(t,usv1_r)
plt.ylabel("r (rad/s)",font,labelpad=18)

plt.xlabel("time (s)", font)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.2) 

# cannot use plt.show() before plt.savefig
# otherwise the fig saved will 
    
plt.savefig("NNtest_velocities.svg",dpi = 600 ,format = "svg",bbox_inches = "tight")