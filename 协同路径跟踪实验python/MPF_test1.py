# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 19:27:48 2022

@author: qys
"""

import numpy as np
import matplotlib.pyplot as plt
from usvmodel import USV1
from nncontroller import NNcontroller1

# initialize
ts = 0.02
tfinal = 300
Ns = tfinal/ts

Nxout = 3
Nusv1 = 16

xout = np.zeros((int(Ns),Nxout))
usv1_out = np.zeros((int(Ns),Nusv1))

# path information
def fxy1(x,y,xo,yo,Ro):
    
    f = x**2+y**2-Ro**2
    fx = 2*x
    fy = 2*y
    ffxx = 2
    ffyy = 2
    ffxy = 0
    
    return f,fx,fy,ffxx,ffyy,ffxy

def fxy2(x,y,xo,yo,Ro):
    
    f = 2*x-y
    fx = np.array([2.])
    fy = -np.array([1.])
    ffxx = 0
    ffyy = 0
    ffxy = 0
    
    return f,fx,fy,ffxx,ffyy,ffxy

def TD1(xf,x,Tf,ts):
    xf_dot  = -(xf-x)/Tf
    xf = ts*xf_dot+xf
    return xf,xf_dot

# simulation start
for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        usv1_nu = np.array([[0],[0],[0]]) # usv1-usv3 states initialization
        usv1_eta = np.array([[6],[-3],[0]])
        usv1_Tf = np.array([[0],[0]])
        usv1_Tc = np.array([[50],[50]])
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
        xo = np.array([0.])
        yo = np.array([0.])
        uo = np.array([0.])  # target position initialization
        Ro = np.array([10.])  # desired enclosing radius 
        psio = np.array([0.])
    
    # target
    ro = 0.02*np.sin(0.02*t)
    xo_dot = uo*np.cos(psio)
    yo_dot = uo*np.sin(psio)
    psio_dot = ro
    xo = ts*xo_dot+xo
    yo = ts*yo_dot+yo
    psio = ts*psio_dot+psio
        
    # deisred heading angle
    f,fx,fy,ffxx,ffyy,ffxy = fxy1(usv1_eta[0],usv1_eta[1],xo,yo,Ro)
    usv_x = usv1_eta[0]
    usv_y = usv1_eta[1]
    if fx>=0 and fy<0:
        chid = -np.arctan(fx/fy)
    if fx>0 and fy>0:
        chid = np.pi-np.arctan(fx/fy)
    if fx<=0 and fy>0:
        chid = -np.arctan(fx/fy)-np.pi
    if fx<0 and fy<0:
        chid = -np.arctan(fx/fy)
    if fy==0:
        if fx>0:
            chid = np.pi/2
        else:
            chid = -np.pi/2
    
    fxy_norm = np.linalg.norm(np.array([fx,fy]))
    
    # get usv states and calculate important parameters
    psi = usv1_eta[2]
    v = usv1_nu[1]
    u = usv1_nu[0]
    betad = np.arctan(v/uc)
    Ud = np.sqrt(uc**2+v**2)
    betad_dot = uc*usv1_nu_dot[1]/(Ud**2)
    chi = psi+betad
    chie = psi-chid+betad
    if chie > np.pi:
        chie = chie % (np.pi*2)-np.pi*2
    if chie<-np.pi:
        chie = chie % (np.pi*2)+np.pi*2
    f_dot = fxy_norm*Ud*np.sin(chie)
    chid_dot = ((fx*ffxy-fy*ffxx)*u*np.cos(chi)+\
                  (fx*ffyy-fy*ffxy)*u*np.sin(chi))/(fxy_norm**2)
    
    # guidance law
    theta = uo*np.sin(psio-chid)
    if t==0:
        k1 = 0.001
        k2 = 0.5
    sigma = (k1*f*(theta**2)+np.sqrt((1+(k1*f)**2)*(theta**2)*(uc**2)-theta**4))/(uc**2-theta**2)
    if t == 0:
        Tf = 0.1
        sigmaf = sigma
    sigmaf,sigmaf_dot = TD1(sigmaf,sigma,Tf,ts) 
    s = chie+np.arctan(k1*f) 
    rc = chid_dot-betad_dot-k2*s
    
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
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
        
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)

    usv1_X = np.array([usv1_nu[0],usv1_nu[1],usv1_nu[0]*usv1_nu[2]])

    usv1_w = np.array([usv1_nu[0],usv1_nu[2]])
     
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)
    
    usv1_ew = usv1_w-usv1_wd-usv1_Lambda
    
    usv1_wdf, usv1_wdf_dot = nnc1.TD(usv1_wd,usv1_wdf,ts)
    
    usv1_Fhat, usv1_What = nnc1.estimator(usv1_X,usv1_What,usv1_ew,ts)
    
    usv1_Tc = nnc1.control_law(usv1_ew,usv1_Fhat,usv1_wdf_dot,usv1_Lambda)
       
    # store time series
    xout[i,:] = np.hstack((np.array([t]),xo[0],yo[0]))
    usv1_out[i,:] = np.hstack((usv1_nu[0],usv1_nu[1],usv1_nu[2],\
                               usv1_eta[0],usv1_eta[1],usv1_eta[2],\
                               usv1_Tc[0],usv1_Tc[1],usv1_Tf[0],\
                               usv1_Tf[1],usv1_fu[0],usv1_fr[0],\
                               usv1_Fhat[0],usv1_Fhat[1],chie[0],chid[0]))
    
    
    
t = xout[:,0]
xo = xout[:,1]
yo = xout[:,2]

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
chie = usv1_out[:,14]
chid = usv1_out[:,15]


# plots
font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12,
}

# north-east map
fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(usv1_y,usv1_x,'r-')
plt.plot(yo,xo,'b--')
thetao = np.arange(0,100,0.1)
xd = Ro*np.cos(thetao)
yd = Ro*np.sin(thetao)
plt.plot(yd,xd)
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

fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(t,chid,'r-',label=r"$\chi_{d}$")
plt.xlabel("time (s)")
plt.ylabel(r"$\chi_{d}$")

fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(t,chie,'r-')
plt.xlabel("time (s)")
plt.ylabel(r"$\chi_{e}$")









