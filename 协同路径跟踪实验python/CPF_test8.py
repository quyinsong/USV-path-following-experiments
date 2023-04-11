# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 18:05:25 2022

@author: qys
"""


# 采用动态面设计协同制导率
# 在CPF_test5基础上增加了事件触发采样控制
# 与CPF_test6的区别是采用ESO估计扰动，而不是NN
# NN估计的收敛太依赖跟踪误差，不适用于事件触发控制

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import  ConnectionPatch
from matplotlib.pyplot import MultipleLocator
from usvmodel import USV1
from esocontroller import ESOcontroller1

############################################################## initialize
ts = 0.05
tfinal = 400
Ns = tfinal/ts

Nxout = 4
Nusv1 = 16
Nusv2 = 16
Nusv3 = 16
Nusv4 = 16
Nusv5 = 16
Nerror = 15
xout = np.zeros((int(Ns),Nxout))
usv1_out = np.zeros((int(Ns),Nusv1))
usv2_out = np.zeros((int(Ns),Nusv2))
usv3_out = np.zeros((int(Ns),Nusv3))
usv4_out = np.zeros((int(Ns),Nusv4))
usv5_out = np.zeros((int(Ns),Nusv5))
error_out = np.zeros((int(Ns),Nerror))

####################################################### path information
def fxy1(x,y):

    f = -10*np.sin(np.pi*x/40)+y-0.5*x
    fx = -0.5-1*np.pi*np.cos(np.pi*x/40)/4
    fy = 1
    ffxx = 1*(np.pi**2)*np.sin(np.pi*x/40)/160
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy

def fxy2(x,y):
    
    f = y-x
    fx = -1
    fy = 1
    ffxx = 0
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy

def fxy3(x,y):
    
    f = y-10
    fx = 0
    fy = 1
    ffxx = 0
    ffyy = 0
    ffxy = 0
    return f,fx,fy,ffxx,ffyy,ffxy

# calculate arc error

def fx1_integrated(x):
    return np.sqrt((-1*np.pi*np.cos(np.pi*x/40)/160-1)**2+1)  # integrated function

def fx2_integrated(x):
    return np.sqrt(2)  # integrated function

def fx3_integrated(x):
    return np.sqrt(1)  # integrated function

def arc_length(xi,xj,N):
    '''
    inputs :
        xi,xj : endpoints of integra interval
        N : number of sample points
    outputs:
        Sum : arc length
    '''
    Sum = 0
    h = (xi-xj)/N
    for i in range(N):
        Sum += h*(fx2_integrated(xi+h*i)+fx2_integrated(xi+h*(i+1)))/2
    
    return Sum               # return arc length
    
################################################ Event triggered mechanism
def ETM1(p,ps,triggeredtimex,triggeredtimey,d,time):
    '''
    

    Parameters
    ----------
    p : TYPE
        DESCRIPTION.
    ps : TYPE
        DESCRIPTION.
    triggeredtimex : TYPE
        DESCRIPTION.
    triggeredtimey : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.

    Returns
    -------
    p : TYPE
        DESCRIPTION.
    ps : TYPE
        DESCRIPTION.
    triggeredtimex : TYPE
        DESCRIPTION.
    triggeredtimey : TYPE
        DESCRIPTION.

    '''
    l = len(triggeredtimex)
    for j1 in range(l):
        # communication ETM
        ek1 = p[0,j1]-ps[0,j1]
        ek2 = p[1,j1]-ps[1,j1]
        if np.abs(ek1) >= d[0,j1]:
            ps[0,j1] = p[0,j1]
            triggeredtimex[j1,time] = j1+1
        else:
            ps[0,j1] = ps[0,j1]
        if np.abs(ek2) >= d[1,j1]:
            ps[1,j1] = p[1,j1]
            triggeredtimey[j1,time] = j1+1
        else:
            ps[1,j1] = ps[1,j1]
    return ps,triggeredtimex,triggeredtimey

def ETM2(w,ws,triggeredtimeU,triggeredtimer,d,time):
    '''
    

    Parameters
    ----------
    w : TYPE
        DESCRIPTION.
    ws : TYPE
        DESCRIPTION.
    triggeredtimeU : TYPE
        DESCRIPTION.
    triggeredtimer : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.

    Returns
    -------
    w : TYPE
        DESCRIPTION.
    ws : TYPE
        DESCRIPTION.
    triggeredtimeU : TYPE
        DESCRIPTION.
    triggeredtimer : TYPE
        DESCRIPTION.

    '''
    l = len(triggeredtimeU)
    for j1 in range(l):
        # ESO ETM
        ek1 = w[0,j1]-ws[0,j1]
        ek2 = w[1,j1]-ws[1,j1]
        if np.abs(ek1) >= d[0,j1]:
            ws[0,j1] = w[0,j1]
            triggeredtimeU[j1,time] = j1+1
        else:
            ws[0,j1] = ws[0,j1]
        if np.abs(ek2) >= d[1,j1]:
            ws[1,j1] = w[1,j1]
            triggeredtimer[j1,time] = j1+1
        else:
            ws[1,j1] = ws[1,j1]
    return ws,triggeredtimeU,triggeredtimer

def ETM3(tauc,taus,triggeredtime1,triggeredtime2,triggeredtime3,d,gamma,time):
    '''
    

    Parameters
    ----------
    tauc : TYPE
        DESCRIPTION.
    taus : TYPE
        DESCRIPTION.
    triggeredtime1 : TYPE
        DESCRIPTION.
    triggeredtime2 : TYPE
        DESCRIPTION.
    triggeredtime3 : TYPE
        DESCRIPTION.
    d : TYPE
        DESCRIPTION.
    gamma : TYPE
        DESCRIPTION.
    time : TYPE
        DESCRIPTION.

    Returns
    -------
    taus : TYPE
        DESCRIPTION.
    triggeredtime1 : TYPE
        DESCRIPTION.
    triggeredtime2 : TYPE
        DESCRIPTION.
    triggeredtime3 : TYPE
        DESCRIPTION.

    '''
    l = len(triggeredtime1)
    for j1 in range(l):
        # surge force event trigger and yaw moment event trigger
        ek1 = tauc[0,j1]-taus[0,j1]
        ek2 = tauc[1,j1]-taus[1,j1]
        if np.abs(ek1) >= gamma[0,j1]*np.abs(tauc[0,j1])+d[0,j1] or np.abs(ek2) >= gamma[1,j1]*np.abs(tauc[1,j1])+d[1,j1]:
            taus[0,j1] = tauc[0,j1]
            taus[1,j1] = tauc[1,j1]
            triggeredtime3[j1,time] = j1+1
        else:
            taus[0,j1] = taus[0,j1]
            taus[1,j1] = taus[1,j1]
        if np.abs(ek1) >= gamma[0,j1]*np.abs(tauc[0,j1])+d[0,j1]:
            triggeredtime1[j1,time] = j1+1
        if np.abs(ek2) >= gamma[1,j1]*np.abs(tauc[1,j1])+d[1,j1]:
            triggeredtime2[j1,time] = j1+1
    return taus,triggeredtime1,triggeredtime2,triggeredtime3
                
    
##################################################### simulation start
for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        usv1_nu = np.array([[0],[0],[0]]) # usv1-usv5 states initialization
        usv1_eta = np.array([[0],[-5],[0]])
        usv1_Tf = np.array([[0],[0]])
        usv1_Ts = np.array([[0],[0]])
        usv2_nu = np.array([[0],[0],[0]])
        usv2_eta = np.array([[0],[-20],[0]])
        usv2_Tf = np.array([[0],[0]])
        usv2_Ts = np.array([[0],[0]])
        usv3_nu = np.array([[0],[0],[0]])
        usv3_eta = np.array([[0],[-30],[0]])
        usv3_Tf = np.array([[0],[0]])
        usv3_Ts = np.array([[0],[0]])
        usv4_nu = np.array([[0],[0],[0]])
        usv4_eta = np.array([[20],[50],[0]])
        usv4_Tf = np.array([[0],[0]])
        usv4_Ts = np.array([[0],[0]])
        usv5_nu = np.array([[0],[0],[0]])
        usv5_eta = np.array([[30],[100],[0]])
        usv5_Tf = np.array([[0],[0]])
        usv5_Ts = np.array([[0],[0]])
        x0 = np.array([0])  # virtual target position initialization
        y0 = np.array([0])
        u0 = np.array([1])
        A = np.array([[0,0,0,0,0],
                      [1,0,0,0,0],
                      [0,1,0,0,0],
                      [1,0,0,0,0],
                      [0,0,0,1,0]]) # communication networks
        B = np.diag([1,0,0,0,0])
        L = np.array([[0,0,0,0,0],
                      [-1,1,0,0,0],
                      [0,-1,1,0,0],
                      [-1,0,0,1,0],
                      [0,0,0,-1,1]])

        Df = np.array([[0,  0,  0, 0,  0],
                       [-16, 0,  0,  0,  0],
                       [0, -16,  0,  0,  0],
                       [16,  0,  0,  0, 0],
                       [0,  0,  0,  16,  0]])  # desired lati-error
        Df0 = np.array([[0],[0],[0],[0],[0]])
        Dl = np.array([[0, 0, 0, 0, 0],
                       [-8*np.sqrt(2), 0, 0, 0, 0],
                       [0, -8*np.sqrt(2), 0, 0, 0],
                       [-24*np.sqrt(2), 0, 0, 0, 0],
                       [0, 0, 0, -16*np.sqrt(2), 0]]) # desired logi-error
        Dl0 = np.array([[0],[0],[0],[0],[0]])
        
        k1 = 0.2*np.ones((5,1)) # guidane law parameters
        k2 = 0.6*np.ones((5,1))
        k3 = 0.1*np.ones((5,1))
    
    # change formation
    if t>=200:
        Df = np.array([[0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0]])  # desired lati-error
        Dl = np.array([[0, 0, 0, 0, 0],
                       [-45, 0, 0, 0, 0],
                       [0, -15, 0, 0, 0],
                       [-15, 0, 0, 0, 0],
                       [0, 0, 0, -15, 0]]) # desired logi-error
    
    # disturbancs and currents
    
    tw = np.array([[5*np.sin(0.1*t)*np.cos(0.05*t)+5],
                   [2*np.sin(0.1*t)*np.cos(0.05*t)+2],
                   [2*np.sin(0.1*t)*np.cos(0.05*t)+2]])
    
    '''
    tw = np.array([[0],[0],[0]])
    '''
    vc = np.array([[0],[0],[0]])
    
    # instantiate USV object
    usv1 = USV1()
    usv2 = USV1()
    usv3 = USV1()
    usv4 = USV1()
    usv5 = USV1()
    usv1_nu,usv1_eta,usv1_nu_dot,usv1_eta_dot,usv1_Tf,usv1_fU,usv1_fr = \
        usv1.state_update2(usv1_nu,usv1_eta,usv1_Tf,usv1_Ts,vc,tw,ts)
    usv2_nu,usv2_eta,usv2_nu_dot,usv2_eta_dot,usv2_Tf,usv2_fU,usv2_fr = \
        usv2.state_update2(usv2_nu,usv2_eta,usv2_Tf,usv2_Ts,vc,tw,ts)
    usv3_nu,usv3_eta,usv3_nu_dot,usv3_eta_dot,usv3_Tf,usv3_fU,usv3_fr = \
        usv3.state_update2(usv3_nu,usv3_eta,usv3_Tf,usv3_Ts,vc,tw,ts)
    usv4_nu,usv4_eta,usv4_nu_dot,usv4_eta_dot,usv4_Tf,usv4_fU,usv4_fr = \
        usv4.state_update2(usv4_nu,usv4_eta,usv4_Tf,usv4_Ts,vc,tw,ts)
    usv5_nu,usv5_eta,usv5_nu_dot,usv5_eta_dot,usv5_Tf,usv5_fU,usv5_fr = \
        usv5.state_update1(usv5_nu,usv5_eta,usv5_Tf,usv5_Ts,vc,tw,ts)
    
    # virtual target state update
    f0,fx0,fy0,ffxx0,ffyy0,ffxy0 = fxy2(x0,y0)
    chi0 = -np.arctan(fx0/fy0)
    x0_dot = u0*np.cos(chi0)
    y0_dot = u0*np.sin(chi0)
    x0 = x0_dot*ts+x0
    y0 = y0_dot*ts+y0
    
    # get all the states of USVs
    eta = np.hstack((usv1_eta,usv2_eta,usv3_eta,usv4_eta,usv5_eta))  # position vector compact form
    nu = np.hstack((usv1_nu,usv2_nu,usv3_nu,usv4_nu,usv5_nu))      # velocity vector compact form
    
    # ETM for communication
    if t==0:
        usv_ps = np.zeros((2,5))
        triggeredtimex = np.zeros((5,int(Ns)))
        triggeredtimey = np.zeros((5,int(Ns)))
        d_xy = np.array([[0.005,0.005,0.005,0.005,0.005],[0.005,0.005,0.005,0.005,0.005]])
    usv_p = eta[0:2,0:5]
    ps,triggeredtimex,triggeredtimey, = ETM1(usv_p,usv_ps,triggeredtimex,triggeredtimey,d_xy,i)
    
    # calculate lati-coordinated errors
    usv_fx = np.zeros((5,1))
    usv_fy = np.zeros((5,1))
    usv_ffxx = np.zeros((5,1))
    usv_ffyy = np.zeros((5,1))
    usv_ffxy = np.zeros((5,1))
    
    F = np.zeros((5,1))   # request memory to store lati-coordinated errors
    
    for j in range(5):
        F[j],usv_fx[j],usv_fy[j],usv_ffxx[j],usv_ffyy[j],usv_ffxy[j] = fxy2(usv_p[0,j],usv_p[1,j])
    
    sum1 = np.sum(A*Df,axis=1)
    Ef = np.dot((L+B),F)-np.array([[sum1[0]],[sum1[1]],[sum1[2]],[sum1[3]],[sum1[4]]])-np.dot(B,Df0)
            
    # calculate longi-coordinated errors
    Lij = np.zeros((5,5))   # request memory to store logi-coordinated errors
    Li0 = np.zeros((5,1))
    
    for m in range(5):
        for n in range(5):
            Lij[m,n] = arc_length(usv_p[0,m],usv_p[0,n],50)
            Li0[m,0] = arc_length(usv_p[0,m],x0,50)
    
    sum2 = np.sum(A*(Lij-Dl),axis=1)      
    El = np.array([[sum2[0]],[sum2[1]],[sum2[2]],[sum2[3]],[sum2[4]]])+np.dot(B,(Li0-Dl0))
    
    # calculate essential parameters
    chid = np.zeros((5,1))
    chid_dot = np.zeros((5,1))
    fxynorm2 = np.zeros((5,1))
    betad = np.zeros((5,1))
    chie = np.zeros((5,1))
    F_dot = np.zeros((5,1))
    
    usv_U = np.array(([np.sqrt(usv1_nu[0]**2+usv1_nu[1]**2),
                       np.sqrt(usv2_nu[0]**2+usv2_nu[1]**2),
                       np.sqrt(usv3_nu[0]**2+usv3_nu[1]**2),
                       np.sqrt(usv4_nu[0]**2+usv4_nu[1]**2),
                       np.sqrt(usv5_nu[0]**2+usv5_nu[1]**2)]))
    usv_Ud = np.array([np.sqrt(u0**2+usv1_nu[1]**2),
                       np.sqrt(u0**2+usv2_nu[1]**2),
                       np.sqrt(u0**2+usv3_nu[1]**2),
                       np.sqrt(u0**2+usv4_nu[1]**2),
                       np.sqrt(u0**2+usv5_nu[1]**2)])
    usv_v = np.array([[nu[1,0]],[nu[1,1]],[nu[1,2]],[nu[1,3]],[nu[1,4]]])
    usv_betad_dot = np.array([u0*usv1_nu_dot[1],
                              u0*usv2_nu_dot[1],
                              u0*usv3_nu_dot[1],
                              u0*usv4_nu_dot[1],
                              u0*usv5_nu_dot[1]])/(np.sqrt(u0**2-usv_v**2))
    
    for j in range(5):
        chid[j] = -np.arctan(usv_fx[j]/usv_fy[j])
        fxynorm2[j] = np.square(usv_fx[j])+np.square(usv_fy[j])
        chid_dot[j] = ((usv_fx[j]*usv_ffxy[j]-usv_fy[j]*usv_ffxx[j])*nu[0,j]*np.cos(eta[2,j])+
                       (usv_fx[j]*usv_ffyy[j]-usv_fy[j]*usv_ffxy[j])*nu[0,j]*np.sin(eta[2,j]))/fxynorm2[j]
        betad[j] = np.arcsin(nu[1,j]/u0)
        chie[j] = eta[2,j]+betad[j]-chid[j]
        F_dot[j] = np.sqrt(fxynorm2[j])*(nu[0,j]*np.sin(eta[2,j]-chid[j])+
                                         nu[1,j]*np.cos(eta[2,j]-chid[j]))
    
    Ef_dot = np.dot((L+B),F_dot)
    
    # coordinated guidance law
    if t == 0:
        psicf = np.zeros((5,1))
        psicf_dot = np.zeros((5,1))
        Uc = np.zeros((5,1))
        rc = np.zeros((5,1))
        psic = np.zeros((5,1))
        psic_dot = np.zeros((5,1))
        e_psi = np.zeros((5,1))
    
    if t == 0:
        Tf = 0.5
    
    for j in range(5):
        psic[j] = chid[j]-betad[j]-np.arctan(k1[j]*Ef[j])
        psicf_dot[j] = -(psicf[j]-psic[j])/Tf
        psicf[j] = ts*psicf_dot[j]+psicf[j]
        psic_dot[j] = chid_dot[j]-usv_betad_dot[j]-\
                      k1[j]*(Ef_dot[j])/(1+np.square(k1[j]*Ef[j]))
        Uc[j] = u0*(1-np.tanh(k3[j]*El[j]))
        e_psi[j] = eta[2,j]-psic[j]
        if e_psi[j]>=np.pi:
            e_psi[j] = e_psi[j] % (np.pi*2)-np.pi*2
        if e_psi[j]<=-np.pi:
            e_psi[j] = e_psi[j] % (np.pi*2)+np.pi*2
        rc[j] = psicf_dot[j]-k2[j]*e_psi[j]
    
    # ESO-based controller
    usv1_w = np.array([usv_U[0],usv1_nu[2]])
    usv2_w = np.array([usv_U[1],usv2_nu[2]])
    usv3_w = np.array([usv_U[2],usv3_nu[2]])
    usv4_w = np.array([usv_U[3],usv4_nu[2]])
    usv5_w = np.array([usv_U[4],usv5_nu[2]])  # w
    
    usv1_wd = np.array([Uc[0],rc[0]])
    usv2_wd = np.array([Uc[1],rc[1]])
    usv3_wd = np.array([Uc[2],rc[2]])
    usv4_wd = np.array([Uc[3],rc[3]])
    usv5_wd = np.array([Uc[4],rc[4]])
    
    usv1_Delta_T = usv1_Ts-usv1_Tf
    usv2_Delta_T = usv2_Ts-usv2_Tf
    usv3_Delta_T = usv3_Ts-usv3_Tf
    usv4_Delta_T = usv4_Ts-usv4_Tf
    usv5_Delta_T = usv5_Ts-usv5_Tf
    
    if t==0:
        B_usv = np.array([[1/50.5,1/50.5],[0.26/17.21,-0.26/17.21]])
        B1_usv = np.array([[1,1],[0.26,-0.26]])
        M_usv = np.array([[50.5,0],[0,17.21]])
        usv1_wdf = usv1_wd
        usv2_wdf = usv2_wd
        usv3_wdf = usv3_wd
        usv4_wdf = usv4_wd
        usv5_wdf = usv5_wd
        usv1_Lambda = np.dot(B_usv,usv1_Delta_T)
        usv2_Lambda = np.dot(B_usv,usv2_Delta_T)
        usv3_Lambda = np.dot(B_usv,usv3_Delta_T)
        usv4_Lambda = np.dot(B_usv,usv4_Delta_T)
        usv5_Lambda = np.dot(B_usv,usv5_Delta_T)
        usv1_what = usv1_w
        usv2_what = usv2_w
        usv3_what = usv3_w
        usv4_what = usv4_w
        usv5_what = usv5_w
        usv1_Fhat = np.array([[0.],[0.]])
        usv2_Fhat = np.array([[0.],[0.]])
        usv3_Fhat = np.array([[0.],[0.]])
        usv4_Fhat = np.array([[0.],[0.]])
        usv5_Fhat = np.array([[0.],[0.]])
        
    esoc1 = ESOcontroller1(Ko1=np.diag([2,2]),Ko2=np.diag([2,2]),\
                           KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    esoc2 = ESOcontroller1(Ko1=np.diag([2,2]),Ko2=np.diag([2,2]),\
                           KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    esoc3 = ESOcontroller1(Ko1=np.diag([2,2]),Ko2=np.diag([2,2]),\
                           KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    esoc4 = ESOcontroller1(Ko1=np.diag([2,2]),Ko2=np.diag([2,2]),\
                           KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    esoc5 = ESOcontroller1(Ko1=np.diag([2,2]),Ko2=np.diag([2,2]),\
                           KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    
    # auxiliary system
    usv1_Lambda = esoc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)
    usv2_Lambda = esoc2.auxiliary_system(usv2_Lambda,usv2_Delta_T,ts)
    usv3_Lambda = esoc3.auxiliary_system(usv3_Lambda,usv3_Delta_T,ts)
    usv4_Lambda = esoc4.auxiliary_system(usv4_Lambda,usv4_Delta_T,ts)
    usv5_Lambda = esoc5.auxiliary_system(usv5_Lambda,usv5_Delta_T,ts)
    
    
    usv1_ew = usv1_w-usv1_wd-usv1_Lambda
    usv2_ew = usv2_w-usv2_wd-usv2_Lambda
    usv3_ew = usv3_w-usv3_wd-usv3_Lambda
    usv4_ew = usv4_w-usv4_wd-usv4_Lambda
    usv5_ew = usv5_w-usv5_wd-usv5_Lambda  # tracking errors
    
    usv1_wdf, usv1_wdf_dot = esoc1.TD(usv1_wd,usv1_wdf,ts)
    usv2_wdf, usv2_wdf_dot = esoc2.TD(usv2_wd,usv2_wdf,ts)
    usv3_wdf, usv3_wdf_dot = esoc3.TD(usv3_wd,usv3_wdf,ts)
    usv4_wdf, usv4_wdf_dot = esoc4.TD(usv4_wd,usv4_wdf,ts)
    usv5_wdf, usv5_wdf_dot = esoc5.TD(usv5_wd,usv5_wdf,ts)
    
    # ETM for ET-ESO
    if t==0:
        usv_ws = np.zeros((2,5))
        triggeredtimeU = np.zeros((5,int(Ns)))
        triggeredtimer = np.zeros((5,int(Ns)))
        d_w = np.array([[0.01,0.01,0.01,0.01,0.01],[0.01,0.01,0.01,0.01,0.01]])
    usv_w = np.hstack((usv1_w,usv2_w,usv3_w,usv4_w,usv5_w))
    
    usv_ws,triggeredtimeU,triggeredtimer, = ETM2(usv_w,usv_ws,triggeredtimeU,triggeredtimer,d_w,i)
    
    usv1_ws = np.array([usv_ws[0:2,0]]).T
    usv2_ws = np.array([usv_ws[0:2,1]]).T
    usv3_ws = np.array([usv_ws[0:2,2]]).T
    usv4_ws = np.array([usv_ws[0:2,3]]).T
    usv5_ws = np.array([usv_ws[0:2,4]]).T
    
    # estimate disturbances by ET-ESO
    usv1_Fhat, usv1_what = esoc1.estimator(usv1_what,usv1_ws,usv1_Fhat,usv1_Tf,ts)
    usv2_Fhat, usv2_what = esoc2.estimator(usv2_what,usv2_ws,usv2_Fhat,usv2_Tf,ts)
    usv3_Fhat, usv3_what = esoc3.estimator(usv3_what,usv3_ws,usv3_Fhat,usv3_Tf,ts)
    usv4_Fhat, usv4_what = esoc4.estimator(usv4_what,usv4_ws,usv4_Fhat,usv4_Tf,ts)
    usv5_Fhat, usv5_what = esoc5.estimator(usv5_what,usv5_ws,usv5_Fhat,usv5_Tf,ts)
    
    # real disturbances
    usv1_F = np.array([usv1_fU,usv1_fr])
    usv2_F = np.array([usv2_fU,usv2_fr])
    usv3_F = np.array([usv3_fU,usv3_fr])
    usv4_F = np.array([usv4_fU,usv4_fr])
    usv5_F = np.array([usv5_fU,usv5_fr])

    # controller
    usv1_Tc = esoc1.control_law(usv1_ew,usv1_Fhat,usv1_wdf_dot,usv1_Lambda)
    usv2_Tc = esoc2.control_law(usv2_ew,usv2_Fhat,usv2_wdf_dot,usv2_Lambda)
    usv3_Tc = esoc3.control_law(usv3_ew,usv3_Fhat,usv3_wdf_dot,usv3_Lambda)
    usv4_Tc = esoc4.control_law(usv4_ew,usv4_Fhat,usv4_wdf_dot,usv4_Lambda)
    usv5_Tc = esoc5.control_law(usv5_ew,usv5_Fhat,usv5_wdf_dot,usv5_Lambda)
    
    # ETM for controller
    Tc = np.hstack((usv1_Tc,usv2_Tc,usv3_Tc,usv4_Tc,usv5_Tc))
        
    tauc = np.dot(B1_usv,Tc)
    
    usv1_tauc = np.array([tauc[:,0]]).T
    usv2_tauc = np.array([tauc[:,1]]).T
    usv3_tauc = np.array([tauc[:,2]]).T
    usv4_tauc = np.array([tauc[:,3]]).T
    usv5_tauc = np.array([tauc[:,4]]).T
    
    if t == 0:
        taus = np.zeros((2,5))
        triggeredtime1 = np.zeros((5,int(Ns)))
        triggeredtime2 = np.zeros((5,int(Ns)))
        triggeredtime3 = np.zeros((5,int(Ns)))
        d_coef = np.array([[0.5],[0.5]])
        d     = d_coef*np.ones((2,5))
        gamma_coef = np.array([[0.01],[0.01]])
        gamma = gamma_coef*np.ones((2,5))
        
    tauc1 = (1+0)*tauc
    Tc1 = np.dot(np.linalg.inv(B1_usv),tauc1)
           
    taus,triggeredtime1,triggeredtime2,triggeredtime3 = ETM3(tauc1,taus,triggeredtime1,triggeredtime2,triggeredtime3,d,gamma,time=i)
    
    Ts = np.dot(np.linalg.inv(B1_usv),taus)
    
    usv1_Ts = np.array([Ts[0:2,0]]).T
    usv2_Ts = np.array([Ts[0:2,1]]).T
    usv3_Ts = np.array([Ts[0:2,2]]).T
    usv4_Ts = np.array([Ts[0:2,3]]).T
    usv5_Ts = np.array([Ts[0:2,4]]).T
    
    usv1_taus = np.array([taus[:,0]]).T
    usv2_taus = np.array([taus[:,1]]).T
    usv3_taus = np.array([taus[:,2]]).T
    usv4_taus = np.array([taus[:,3]]).T
    usv5_taus = np.array([taus[:,4]]).T
    
    # store time series
    xout[i,:] = np.hstack((np.array([t]),x0[0],y0[0],u0[0]))
    usv1_out[i,:] = np.hstack((usv1_nu[0],usv1_nu[1],usv1_nu[2],\
                               usv1_eta[0],usv1_eta[1],usv1_eta[2],\
                               Tc1[0,0],Tc1[1,0],usv1_Ts[0],\
                               usv1_Ts[1],usv1_Tf[0],usv1_Tf[1],\
                               usv1_fU[0],usv1_fr[0],usv1_Fhat[0],\
                               usv1_Fhat[1]))
    usv2_out[i,:] = np.hstack((usv2_nu[0],usv2_nu[1],usv2_nu[2],\
                               usv2_eta[0],usv2_eta[1],usv2_eta[2],\
                               Tc1[0,1],Tc1[1,1],usv2_Ts[0],\
                               usv2_Ts[1],usv2_Tf[0],usv2_Tf[1],\
                               usv2_fU[0],usv2_fr[0],usv2_Fhat[0],\
                               usv2_Fhat[1]))
    usv3_out[i,:] = np.hstack((usv3_nu[0],usv3_nu[1],usv3_nu[2],\
                               usv3_eta[0],usv3_eta[1],usv3_eta[2],\
                               Tc1[0,2],Tc1[1,2],usv3_Ts[0],\
                               usv3_Ts[1],usv3_Tf[0],usv3_Tf[1],\
                               usv3_fU[0],usv3_fr[0],usv3_Fhat[0],\
                               usv3_Fhat[1]))
    usv4_out[i,:] = np.hstack((usv4_nu[0],usv4_nu[1],usv4_nu[2],\
                               usv4_eta[0],usv4_eta[1],usv4_eta[2],\
                               Tc1[0,3],Tc1[1,3],usv4_Ts[0],\
                               usv4_Ts[1],usv4_Tf[0],usv4_Tf[1],\
                               usv4_fU[0],usv4_fr[0],usv4_Fhat[0],\
                               usv4_Fhat[1]))
    usv5_out[i,:] = np.hstack((usv5_nu[0],usv5_nu[1],usv5_nu[2],\
                               usv5_eta[0],usv5_eta[1],usv5_eta[2],\
                               Tc1[0,4],Tc1[1,4],usv5_Ts[0],\
                               usv5_Ts[1],usv5_Tf[0],usv5_Tf[1],\
                               usv5_fU[0],usv5_fr[0],usv5_Fhat[0],\
                               usv5_Fhat[1]))
    error_out[i,:] = np.hstack((Ef[0],Ef[1],Ef[2],Ef[3],Ef[4],\
                                El[0],El[1],El[2],El[3],El[4],\
                                chie[0],chie[1],chie[2],chie[3],chie[4]))
    
    
################################################################ get data 
t = xout[:,0]
x0 = xout[:,1]
y0 = xout[:,2]
u0 = xout[:,3]

usv1_u = usv1_out[:,0]
usv1_v = usv1_out[:,1]
usv1_r = usv1_out[:,2]
usv1_x = usv1_out[:,3]
usv1_y = usv1_out[:,4]
usv1_psi = usv1_out[:,5]
usv1_Tc0_1 = usv1_out[:,6]
usv1_Tc0_2 = usv1_out[:,7]
usv1_Tc1_1 = usv1_out[:,8]
usv1_Tc1_2 = usv1_out[:,9]
usv1_Tf1 = usv1_out[:,10]
usv1_Tf2 = usv1_out[:,11]
usv1_fU = usv1_out[:,12]
usv1_fr = usv1_out[:,13]
usv1_fUhat = usv1_out[:,14]
usv1_frhat = usv1_out[:,15]
usv1_U = np.sqrt(usv1_u**2+usv1_v**2)

usv2_u = usv2_out[:,0]
usv2_v = usv2_out[:,1]
usv2_r = usv2_out[:,2]
usv2_x = usv2_out[:,3]
usv2_y = usv2_out[:,4]
usv2_psi = usv2_out[:,5]
usv2_Tc0_1 = usv2_out[:,6]
usv2_Tc0_2 = usv2_out[:,7]
usv2_Tc1_1 = usv2_out[:,8]
usv2_Tc1_2 = usv2_out[:,9]
usv2_Tf1 = usv2_out[:,10]
usv2_Tf2 = usv2_out[:,11]
usv2_fU = usv2_out[:,12]
usv2_fr = usv2_out[:,13]
usv2_fUhat = usv2_out[:,14]
usv2_frhat = usv2_out[:,15]
usv2_U = np.sqrt(usv2_u**2+usv2_v**2)

usv3_u = usv3_out[:,0]
usv3_v = usv3_out[:,1]
usv3_r = usv3_out[:,2]
usv3_x = usv3_out[:,3]
usv3_y = usv3_out[:,4]
usv3_psi = usv3_out[:,5]
usv3_Tc0_1 = usv3_out[:,6]
usv3_Tc0_2 = usv3_out[:,7]
usv3_Tc1_1 = usv3_out[:,8]
usv3_Tc1_2 = usv3_out[:,9]
usv3_Tf1 = usv3_out[:,10]
usv3_Tf2 = usv3_out[:,11]
usv3_fU = usv3_out[:,12]
usv3_fr = usv3_out[:,13]
usv3_fUhat = usv3_out[:,14]
usv3_frhat = usv3_out[:,15]
usv3_U = np.sqrt(usv3_u**2+usv3_v**2)

usv4_u = usv4_out[:,0]
usv4_v = usv4_out[:,1]
usv4_r = usv4_out[:,2]
usv4_x = usv4_out[:,3]
usv4_y = usv4_out[:,4]
usv4_psi = usv4_out[:,5]
usv4_Tc0_1 = usv4_out[:,6]
usv4_Tc0_2 = usv4_out[:,7]
usv4_Tc1_1 = usv4_out[:,8]
usv4_Tc1_2 = usv4_out[:,9]
usv4_Tf1 = usv4_out[:,10]
usv4_Tf2 = usv4_out[:,11]
usv4_fu = usv4_out[:,12]
usv4_fr = usv4_out[:,13]
usv4_fUhat = usv4_out[:,14]
usv4_frhat = usv4_out[:,15]
usv4_U = np.sqrt(usv4_u**2+usv4_v**2)

usv5_u = usv5_out[:,0]
usv5_v = usv5_out[:,1]
usv5_r = usv5_out[:,2]
usv5_x = usv5_out[:,3]
usv5_y = usv5_out[:,4]
usv5_psi = usv5_out[:,5]
usv5_Tc0_1 = usv5_out[:,6]
usv5_Tc0_2 = usv5_out[:,7]
usv5_Tc1_1 = usv5_out[:,8]
usv5_Tc1_2 = usv5_out[:,9]
usv5_Tf1 = usv5_out[:,8]
usv5_Tf2 = usv5_out[:,9]
usv5_fU = usv5_out[:,10]
usv5_fr = usv5_out[:,11]
usv5_fUhat = usv5_out[:,12]
usv5_frhat = usv5_out[:,13]
usv5_U = np.sqrt(usv5_u**2+usv5_v**2)

ef1 = error_out[:,0]
ef2 = error_out[:,1]
ef3 = error_out[:,2]
ef4 = error_out[:,3]
ef5 = error_out[:,4]
el1 = error_out[:,5]
el2 = error_out[:,6]
el3 = error_out[:,7]
el4 = error_out[:,8]
el5 = error_out[:,9]
chie1 = error_out[:,10]
chie2 = error_out[:,11]
chie3 = error_out[:,12]
chie4 = error_out[:,13]
chie5 = error_out[:,14]

############################################# plots

font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12,
}   # font set

linewid = 1
################ north-east map
def vhicleplot(usv_len,pos,psi,linest,linewid):
    Rbn = np.array([[np.cos(psi),-np.sin(psi)],
                    [np.sin(psi),np.cos(psi)]])
    '''
       C  
     B   D
     A   E
    '''
    length = usv_len
    width = length*2/5
    Ab = np.array([[-length*3/5],[-width/2]])
    Bb = np.array([[0.],[-width/2]])
    Cb = np.array([[length*2/5],[0]])
    Db = np.array([[0.],[width/2]])
    Eb = np.array([[-length*3/5],[width/2]])
    
    AI = np.dot(Rbn,Ab)+pos
    BI = np.dot(Rbn,Bb)+pos
    CI = np.dot(Rbn,Cb)+pos
    DI = np.dot(Rbn,Db)+pos
    EI = np.dot(Rbn,Eb)+pos
    
    xI = np.array([AI[0],BI[0],CI[0],DI[0],EI[0],AI[0]])
    yI = np.array([AI[1],BI[1],CI[1],DI[1],EI[1],AI[1]])
    plt.plot(yI,xI,linest,linewidth=linewid)
    
fig = plt.figure(figsize=(6,6),dpi = 600)
plt.plot(usv1_y,usv1_x,'r-',linewidth = linewid,label="vehicle 1")
plt.plot(usv2_y,usv2_x,'g-',linewidth = linewid,label="vehicle 2")
plt.plot(usv3_y,usv3_x,'b-',linewidth = linewid,label="vehicle 3")
plt.plot(usv4_y,usv4_x,'c-',linewidth = linewid,label="vehicle 4")
plt.plot(usv5_y,usv5_x,'m-',linewidth = linewid,label="vehicle 5")
plt.plot(y0,x0,'k--',linewidth = linewid)

usv_len = 5.
for j in np.arange(0,int(Ns),2000):
    vhicleplot(usv_len=usv_len,pos=np.array([[usv1_x[j]],[usv1_y[j]]]),\
               psi=float(usv1_psi[j]),linest = 'r-',linewid = linewid)
    vhicleplot(usv_len,pos=np.array([[usv2_x[j]],[usv2_y[j]]]),\
               psi = usv2_psi[j],linest = 'g-',linewid = linewid)
    vhicleplot(usv_len,pos=np.array([[usv3_x[j]],[usv3_y[j]]]),\
               psi = usv3_psi[j],linest = 'b-',linewid = linewid)
    vhicleplot(usv_len,pos=np.array([[usv4_x[j]],[usv4_y[j]]]),\
               psi = usv4_psi[j],linest = 'c-',linewid = linewid)
    vhicleplot(usv_len,pos=np.array([[usv5_x[j]],[usv5_y[j]]]),\
               psi = usv5_psi[j],linest = 'm-',linewid = linewid)

plt.legend(loc="right")
plt.xlabel("east (m)",font,labelpad=8.5)
plt.ylabel("north (m)",font,labelpad=8.5)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2map.png",dpi=600,bbox_inches = "tight")

#################### thrust
def zone_and_linked(ax,axins,zone_left,zone_right,x,y,linked='bottom',
                    x_ratio=0.05,y_ratio=0.05):
    """缩放内嵌图形，并且进行连线
    ax:         调用plt.subplots返回的画布。例如： fig,ax = plt.subplots(1,1)
    axins:      内嵌图的画布。 例如 axins = ax.inset_axes((0.4,0.1,0.4,0.3))
    zone_left:  要放大区域的横坐标左端点
    zone_right: 要放大区域的横坐标右端点
    x:          X轴标签
    y:          列表，所有y值
    linked:     进行连线的位置，{'bottom','top','left','right'}
    x_ratio:    X轴缩放比例
    y_ratio:    Y轴缩放比例
    """
    xlim_left = x[zone_left]-(x[zone_right]-x[zone_left])*x_ratio
    xlim_right = x[zone_right]+(x[zone_right]-x[zone_left])*x_ratio

    y_data = np.hstack([yi[zone_left:zone_right] for yi in y])
    ylim_bottom = np.min(y_data)-(np.max(y_data)-np.min(y_data))*y_ratio
    ylim_top = np.max(y_data)+(np.max(y_data)-np.min(y_data))*y_ratio

    axins.set_xlim(xlim_left, xlim_right)
    axins.set_ylim(ylim_bottom, ylim_top)

    ax.plot([xlim_left,xlim_right,xlim_right,xlim_left,xlim_left],
            [ylim_bottom,ylim_bottom,ylim_top,ylim_top,ylim_bottom],"black")

    if linked == 'bottom':
        xyA_1, xyB_1 = (xlim_left,ylim_top), (xlim_left,ylim_bottom)
        xyA_2, xyB_2 = (xlim_right,ylim_top), (xlim_right,ylim_bottom)
    elif  linked == 'top':
        xyA_1, xyB_1 = (xlim_left,ylim_bottom), (xlim_left,ylim_top)
        xyA_2, xyB_2 = (xlim_right,ylim_bottom), (xlim_right,ylim_top)
    elif  linked == 'left':
        xyA_1, xyB_1 = (xlim_right,ylim_top), (xlim_left,ylim_top)
        xyA_2, xyB_2 = (xlim_right,ylim_bottom), (xlim_left,ylim_bottom)
    elif  linked == 'right':
        xyA_1, xyB_1 = (xlim_left,ylim_top), (xlim_right,ylim_top)
        xyA_2, xyB_2 = (xlim_left,ylim_bottom), (xlim_right,ylim_bottom)
        
    con = ConnectionPatch(xyA=xyA_1,xyB=xyB_1,coordsA="data",
                          coordsB="data",axesA=axins,axesB=ax)
    axins.add_artist(con)
    con = ConnectionPatch(xyA=xyA_2,xyB=xyB_2,coordsA="data",
                          coordsB="data",axesA=axins,axesB=ax)
    axins.add_artist(con)

fig = plt.figure(figsize=(6,6),dpi = 600)
f, ax = plt.subplots(2,1,sharex = True)

plt.subplot(211)
ax[0].plot(t,usv1_Tc0_1,'r-',linewidth = linewid,label=r"$T_{cl1}$")
ax[0].plot(t,usv1_Tc1_1,'g-',linewidth = linewid,label=r"$T_{sl1}$")
plt.legend(loc="upper right")
plt.ylim(-80,180)
plt.xlim(0,tfinal)
plt.ylabel("thrust of vehicle 1 (N)")
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

# 绘制缩放图
axins = ax[0].inset_axes((0.4, 0.2, 0.4, 0.3))
axins.plot(t,usv1_Tc0_1,'r-',linewidth = linewid)
axins.plot(t,usv1_Tc1_1,'g-',linewidth = linewid)

zone_and_linked(ax[0],axins,zone_left=int(30/0.05),zone_right=int(40/0.05),\
                x=t,y=[usv1_Tc0_1,usv1_Tc1_1],linked='right',\
                x_ratio=0.05,y_ratio=0.05)
'''
axins.set_xlim(20, 30)
axins.set_ylim(70, 100)
'''

plt.subplot(212)
ax[1].plot(t,usv1_Tc0_2,'b-',linewidth = linewid,label=r"$T_{cr1}$")
ax[1].plot(t,usv1_Tc1_2,'m-',linewidth = linewid,label=r"$T_{sr1}$")
plt.legend(loc="upper right")
plt.ylim(-80,180)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
plt.ylabel("thrust of vehicle 1 (N)")
plt.xlabel("time (s)")

# 绘制缩放图
axins = ax[1].inset_axes((0.4, 0.2, 0.4, 0.3))
axins.plot(t,usv1_Tc0_2,'b-',linewidth = linewid,label="Tc0_1")
axins.plot(t,usv1_Tc1_2,'m-',linewidth = linewid,label="Tc1_1")
zone_and_linked(ax[1],axins,zone_left=int(30/0.05),zone_right=int(40/0.05),\
                x=t,y=[usv1_Tc0_2,usv1_Tc1_2],linked='right',\
                x_ratio=0.05,y_ratio=0.05)
'''
axins.set_xlim(25, 30)
axins.set_ylim(70, 100)
'''

'''
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

plt.subplot(511)
axarr[0].plot(t,usv1_Tc0_1,'r-',label="Tc0_1")
axarr[0].plot(t,usv1_Tc0_2,'b-',label="Tc1_2")
axarr[0].plot(t,usv1_Tc1_1,'g-',label="Tc1_1")
axarr[0].plot(t,usv1_Tc1_2,'m-',label="Tc1_2")
plt.title("usv1")
plt.legend(loc="upper right")
plt.ylim(-100,200)

plt.subplot(512)
axarr[1].plot(t,usv2_Tc0_1,'r-')
axarr[1].plot(t,usv2_Tc0_2,'b-')
axarr[1].plot(t,usv2_Tc1_1,'g-')
axarr[1].plot(t,usv2_Tc1_2,'m-')
plt.title("usv2")
plt.ylabel("thrust (N)")
plt.ylim(-100,200)

plt.subplot(513)
axarr[2].plot(t,usv3_Tc0_1,'r-')
axarr[2].plot(t,usv3_Tc0_2,'b-')
axarr[2].plot(t,usv3_Tc1_1,'g-')
axarr[2].plot(t,usv3_Tc1_2,'m-')
plt.title("usv3")
plt.ylim(-100,200)

plt.subplot(514)
axarr[3].plot(t,usv4_Tc0_1,'r-')
axarr[3].plot(t,usv4_Tc0_2,'b-')
axarr[3].plot(t,usv4_Tc1_1,'g-')
axarr[3].plot(t,usv4_Tc1_2,'m-')
plt.title("usv4")
plt.ylim(-100,200)

plt.subplot(515)
axarr[4].plot(t,usv5_Tc0_1,'r-')
axarr[4].plot(t,usv5_Tc0_2,'b-')
axarr[4].plot(t,usv5_Tc1_1,'g-')
axarr[4].plot(t,usv5_Tc1_2,'m-')

plt.title("usv5")
plt.ylabel("thrust (N)")
plt.xlabel("time (s)")
plt.xlim(0,tfinal)
plt.ylim(-100,200)
f.subplots_adjust(hspace=0.8)
'''
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2thrust.png",dpi=600,bbox_inches = "tight")

# velocies
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,usv1_U,'r-',linewidth = linewid,label = "vehicle 1")
axarr[0].plot(t,usv2_U,'g-',linewidth = linewid,label = "vehicle 2")
axarr[0].plot(t,usv3_U,'b-',linewidth = linewid,label = "vehicle 3")
axarr[0].plot(t,usv4_U,'c-',linewidth = linewid,label = "vehicle 4")
axarr[0].plot(t,usv5_U,'m-',linewidth = linewid,label = "vehicle 5")
# axarr[0].plot(t,u0,'k--',label = "U0",linewidth = linewid)
plt.legend(loc="upper right",ncol=1)
plt.ylabel(r"$U_{i}$ (m/s)",font,labelpad=8.5)
y_major_locator=MultipleLocator(0.5)#以每15显示

ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
ax1.yaxis.set_major_locator(y_major_locator)

plt.subplot(312)
axarr[1].plot(t,usv1_v,'r-',linewidth = linewid)
axarr[1].plot(t,usv2_v,'g-',linewidth = linewid)
axarr[1].plot(t,usv3_v,'b-',linewidth = linewid)
axarr[1].plot(t,usv4_v,'c-',linewidth = linewid)
axarr[1].plot(t,usv5_v,'m-',linewidth = linewid)

plt.ylabel(r"$v_{i}$ (m/s)",font,labelpad=8.5)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.subplot(313)
axarr[2].plot(t,usv1_r,'r-',linewidth = linewid)
axarr[2].plot(t,usv2_r,'g-',linewidth = linewid)
axarr[2].plot(t,usv3_r,'b-',linewidth = linewid)
axarr[2].plot(t,usv4_r,'c-',linewidth = linewid)
axarr[2].plot(t,usv5_r,'m-',linewidth = linewid)

plt.ylabel(r"$r_{i}$ (rad/s)",font,labelpad=8.5)

plt.xlabel("time (s)", font,labelpad=8.5)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.5)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2velocities.png",dpi=600,bbox_inches = "tight")

# coordinated following errors
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(2,1,sharex = True)

plt.subplot(211)
axarr[0].plot(t,ef1,'r-',linewidth = linewid,label="$e_{f1}$")
axarr[0].plot(t,ef2,'g-',linewidth = linewid,label="$e_{f2}$")
axarr[0].plot(t,ef3,'b-',linewidth = linewid,label="$e_{f3}$")
axarr[0].plot(t,ef4,'c-',linewidth = linewid,label="$e_{f4}$")
axarr[0].plot(t,ef5,'m-',linewidth = linewid,label="$e_{f5}$")

plt.ylabel(r"$E_{f}$ (m)",font,labelpad=8.5)
plt.legend(loc="upper right",labelspacing=0.4,\
           columnspacing=0.4)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.subplot(212)
axarr[1].plot(t,el1,'r-',linewidth = linewid,label = "$e_{l1}$")
axarr[1].plot(t,el2,'g-',linewidth = linewid,label = "$e_{l2}$")
axarr[1].plot(t,el3,'b-',linewidth = linewid,label = "$e_{l3}$")
axarr[1].plot(t,el4,'c-',linewidth = linewid,label = "$e_{l4}$")
axarr[1].plot(t,el5,'m-',linewidth = linewid,label = "$e_{l5}$")

plt.ylabel(r"$E_{s}$ (m)",font,labelpad=8.5)
plt.legend(loc="upper right",labelspacing=0.4,\
           columnspacing=0.4)
'''
plt.subplot(313)
axarr[2].plot(t,chie1,'r-',linewidth = linewid,label = "$\chi_{e1}$")
axarr[2].plot(t,chie2,'g-',linewidth = linewid,label = "$\chi_{e2}$")
axarr[2].plot(t,chie3,'b-',linewidth = linewid,label = "$\chi_{e3}$")
axarr[2].plot(t,chie4,'c-',linewidth = linewid,label = "$\chi_{e4}$")
axarr[2].plot(t,chie5,'m-',linewidth = linewid,label = "$\chi_{e5}$")

plt.legend(loc="upper right")
plt.ylabel(r"$\theta_{e} (rad)$",font,labelpad=18)
'''
plt.xlabel("time (s)",font,labelpad=8.5)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
f.subplots_adjust(hspace=0.5)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2coordinatederrors.png",dpi=600,bbox_inches = "tight")

# estimation erros

fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(2,1,sharex = True)

plt.subplot(211)
axarr[0].plot(t,usv1_fU,'r-',linewidth = linewid,label = r"$f_{U1}$")
axarr[0].plot(t,usv1_fUhat,'b--',linewidth = linewid,label=r"$\hat{f}_{U1}$")
plt.ylabel(r"$f_{U1}$ and $\hat{f}_{U1}$ (N)",font,labelpad=8.5)
plt.legend(loc="upper right",labelspacing=0.4,\
           columnspacing=0.4)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.subplot(212)
axarr[1].plot(t,usv1_fr,'r-',linewidth = linewid,label=r"$f_{r1}$")
axarr[1].plot(t,usv1_frhat,'b--',linewidth = linewid,label=r"$\hat{f}_{r1}$")
plt.ylabel(r"$f_{r1}$ and $\hat{f}_{r1}$ (Nm)",font,labelpad=8.5)
plt.legend(loc="upper right")
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
y_major_locator=MultipleLocator(0.5)#以每15显示
ax1.yaxis.set_major_locator(y_major_locator)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2estimation.png",dpi=600,bbox_inches = "tight")

'''
# tiggered time for comunication x
fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtimex[j1,j2]:
            plt.plot(t[j2],triggeredtimex[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtimex[j1,j2],line_style[j1],\
             label = title[j1], markersize=0.1)
plt.legend(loc = "upper center",ncol = 5,labelspacing=0.4,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $x_{i}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("triggeredtime0_Us.png",dpi=600,bbox_inches = "tight")

# tiggered time for comunication y
fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtimey[j1,j2]:
            plt.plot(t[j2],triggeredtimey[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtimey[j1,j2],line_style[j1],\
             label = title[j1], markersize=0.1)
plt.legend(loc = "upper center",ncol = 5,labelspacing=0.4,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $y_{i}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("triggeredtime0_Us.png",dpi=600,bbox_inches = "tight")
'''

# triggered time U for ESO
fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtimeU[j1,j2]:
            plt.plot(t[j2],triggeredtimeU[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtimeU[j1,j2],line_style[j1],\
             label = title[j1], markersize=0.1)
plt.legend(loc = "upper center",ncol = 5,labelspacing=0.4,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $U_{i}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)

plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2triggeredtime0_Us.png",dpi=600,bbox_inches = "tight")

# triggered time r for ESO

fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtimer[j1,j2]:
            plt.plot(t[j2],triggeredtimer[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtimer[j1,j2],line_style[j1],\
             label = title[j1], markersize=0.1)
plt.legend(loc = "upper center",ncol=5,labelspacing=0.2,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $r_{i}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax = plt.gca()
ax.tick_params(axis='both',tick1On = True,tick2On=True)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2triggeredtime0_rs.png",dpi=600,bbox_inches = "tight")

# triggered time 1
fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime1[j1,j2]:
            plt.plot(t[j2],triggeredtime1[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtime1[j1,j2],line_style[j1],\
             label = title[j1], markersize=8)
plt.legend(loc = "upper center",ncol=5,labelspacing=0.2,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $\tau_{ui}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2triggeredtime1.png",dpi=600,bbox_inches = "tight")

# triggered time 2
fig = plt.figure(figsize=(6,3),dpi = 600)

label = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime2[j1,j2]:
            plt.plot(t[j2],triggeredtime2[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtime2[j1,j2],line_style[j1],\
             label = label[j1], markersize=8)
plt.legend(loc = "upper center",ncol=5,labelspacing=0.2,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $\tau_{ri}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',tick1On = True,tick2On=True)
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2triggeredtime2.png",dpi=600,bbox_inches = "tight")

# triggered time 3
fig = plt.figure(figsize=(6,3),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r+','gx','b|','c1','m2'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime3[j1,j2]:
            plt.plot(t[j2],triggeredtime3[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtime3[j1,j2],line_style[j1],\
             label = title[j1], markersize=8)
plt.legend(loc = "upper center",ncol=5,labelspacing=0.2,\
           columnspacing=0.4)

plt.xlabel("time (s)",labelpad=8.5)
plt.ylabel(r"Sampling instanes of $T_{ci}$",labelpad=8.5)

plt.ylim(0,7)
plt.xlim(0,tfinal)
ax1 = plt.gca()
ax1.tick_params(axis='both',bottom=True, top=True, left=True, right=True )
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.savefig("Case2triggeredtime3.png",dpi=600,bbox_inches = "tight")

# calculate the total triggered times
sum_triggeredtimes1 = np.zeros((5,1))
sum_triggeredtimes2 = np.zeros((5,1))
sum_triggeredtimes3 = np.zeros((5,1))
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime1[j1,j2]:
            sum_triggeredtimes1[j1,0] +=1
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime2[j1,j2]:
            sum_triggeredtimes2[j1,0] +=1
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime3[j1,j2]:
            sum_triggeredtimes3[j1,0] +=1
           





    
    
    
    
    
    
    
    
    
    
    
    
    
    
