# -*- coding: utf-8 -*-
"""
Created on Sat Oct 29 11:45:01 2022

@author: qys
"""

# 采用动态面设计协同制导率
# 在CPF_test5基础上增加了事件触发采样控制

import numpy as np
import matplotlib.pyplot as plt
from usvmodel import USV1
from nncontroller import NNcontroller1

##############################################################3 initialize
ts = 0.05
tfinal = 200
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

#######################################################33 path information
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
        Sum += h*(fx1_integrated(xi+h*i)+fx1_integrated(xi+h*(i+1)))/2
    
    return Sum               # return arc length
    
################################################ Event triggered mechanism
def ETM1(tauc0,tauc1,triggeredtime1,triggeredtime2,d,gamma,time):
    '''

    Parameters
    ----------
    tauc0 : TYPE
        real time control inputs 2x5
    tauc1 : TYPE
        last time control inputs 2x5
    triggeredtime : TYPE
        record triggered time 5xNs
    d : TYPE
        threshold of triggering 5x1
    gamma : TYPE 5x1
        dynamic threshold coefficient 5x1
    time : TYPE
        current time 5x1

    Returns
    -------
    None.

    '''
    l = len(triggeredtime1)
    for j1 in range(l):
        # surge force event trigger and yaw moment event trigger
        ek1 = tauc0[0,j1]-tauc1[0,j1]
        ek2 = tauc0[1,j1]-tauc1[1,j1]
        if np.abs(ek1) >= gamma[0,j1]*np.abs(tauc1[0,j1])+d[0,j1] or np.abs(ek2) >= gamma[1,j1]*np.abs(tauc1[1,j1])+d[1,j1]:
            tauc1[0,j1] = tauc0[0,j1]
            tauc1[1,j1] = tauc0[1,j1]
            triggeredtime1[j1,time] = j1+1
            triggeredtime2[j1,time] = 1
        else:
            tauc1[0,j1] = tauc1[0,j1]
            tauc1[1,j1] = tauc1[1,j1]
            
    return tauc1,triggeredtime1,triggeredtime2
                
    
##################################################### simulation start
for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        usv1_nu = np.array([[0],[0],[0]]) # usv1-usv5 states initialization
        usv1_eta = np.array([[0],[-5],[0]])
        usv1_Tf = np.array([[0],[0]])
        usv1_Tc1 = np.array([[0],[0]])
        usv1_Tc0 = np.array([[0],[0]])
        usv2_nu = np.array([[0],[0],[0]])
        usv2_eta = np.array([[0],[-10],[0]])
        usv2_Tf = np.array([[0],[0]])
        usv2_Tc1 = np.array([[0],[0]])
        usv2_Tc0 = np.array([[0],[0]])
        usv3_nu = np.array([[0],[0],[0]])
        usv3_eta = np.array([[0],[-20],[0]])
        usv3_Tf = np.array([[0],[0]])
        usv3_Tc1 = np.array([[0],[0]])
        usv3_Tc0 = np.array([[0],[0]])
        usv4_nu = np.array([[0],[0],[0]])
        usv4_eta = np.array([[-5],[15],[0]])
        usv4_Tf = np.array([[0],[0]])
        usv4_Tc1 = np.array([[0],[0]])
        usv4_Tc0 = np.array([[0],[0]])
        usv5_nu = np.array([[0],[0],[0]])
        usv5_eta = np.array([[-10],[20],[0]])
        usv5_Tf = np.array([[0],[0]])
        usv5_Tc1 = np.array([[0],[0]])
        usv5_Tc0 = np.array([[0],[0]])
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

        Df = np.array([[0,  8,  0, -8,  0],
                       [-8, 0,  8,  0,  0],
                       [0, -8,  0,  0,  0],
                       [8,  0,  0,  0, -8],
                       [0,  0,  0,  8,  0]])  # desired lati-error
        Df0 = np.array([[0],[0],[0],[0],[0]])
        Dl = np.array([[0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0],
                       [0, 0, 0, 0, 0]]) # desired logi-error
        Dl0 = np.array([[0],[0],[0],[0],[0]])
        
        k1 = 0.2*np.ones((5,1)) # guidane law parameters
        k2 = 0.2*np.ones((5,1))
        k3 = 0.1*np.ones((5,1))
        
    # disturbancs and currents
    
    tw = np.array([[10*np.sin(0.1*t)*np.cos(0.05*t)+5],
                   [5*np.sin(0.1*t)*np.cos(0.05*t)+5],
                   [5*np.sin(0.1*t)*np.cos(0.05*t)+5]])
    
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
        usv1.state_update2(usv1_nu,usv1_eta,usv1_Tf,usv1_Tc1,vc,tw,ts)
    usv2_nu,usv2_eta,usv2_nu_dot,usv2_eta_dot,usv2_Tf,usv2_fU,usv2_fr = \
        usv2.state_update2(usv2_nu,usv2_eta,usv2_Tf,usv2_Tc1,vc,tw,ts)
    usv3_nu,usv3_eta,usv3_nu_dot,usv3_eta_dot,usv3_Tf,usv3_fU,usv3_fr = \
        usv3.state_update2(usv3_nu,usv3_eta,usv3_Tf,usv3_Tc1,vc,tw,ts)
    usv4_nu,usv4_eta,usv4_nu_dot,usv4_eta_dot,usv4_Tf,usv4_fU,usv4_fr = \
        usv4.state_update2(usv4_nu,usv4_eta,usv4_Tf,usv4_Tc1,vc,tw,ts)
    usv5_nu,usv5_eta,usv5_nu_dot,usv5_eta_dot,usv5_Tf,usv5_fU,usv5_fr = \
        usv5.state_update1(usv5_nu,usv5_eta,usv5_Tf,usv5_Tc1,vc,tw,ts)
    
    # virtual target state update
    f0,fx0,fy0,ffxx0,ffyy0,ffxy0 = fxy1(x0,y0)
    chi0 = -np.arctan(fx0/fy0)
    x0_dot = u0*np.cos(chi0)
    y0_dot = u0*np.sin(chi0)
    x0 = x0_dot*ts+x0
    y0 = y0_dot*ts+y0
    
    # get all the states of USVs
    eta = np.hstack((usv1_eta,usv2_eta,usv3_eta,usv4_eta,usv5_eta))  # position vector compact form
    nu = np.hstack((usv1_nu,usv2_nu,usv3_nu,usv4_nu,usv5_nu))      # velocity vector compact form
    
    # calculate lati-coordinated errors
    usv_fx = np.zeros((5,1))
    usv_fy = np.zeros((5,1))
    usv_ffxx = np.zeros((5,1))
    usv_ffyy = np.zeros((5,1))
    usv_ffxy = np.zeros((5,1))
    
    F = np.zeros((5,1))   # request memory to store lati-coordinated errors
    
    for j in range(5):
        F[j],usv_fx[j],usv_fy[j],usv_ffxx[j],usv_ffyy[j],usv_ffxy[j] = fxy1(eta[0,j],eta[1,j])
    
    sum1 = np.sum(A*Df,axis=1)
    Ef = np.dot((L+B),F)-np.array([[sum1[0]],[sum1[1]],[sum1[2]],[sum1[3]],[sum1[4]]])-np.dot(B,Df0)
            
    # calculate longi-coordinated errors
    Lij = np.zeros((5,5))   # request memory to store logi-coordinated errors
    Li0 = np.zeros((5,1))
    
    for m in range(5):
        for n in range(5):
            Lij[m,n] = arc_length(eta[0,m],eta[0,n],50)
            Li0[m,0] = arc_length(eta[0,m],x0,50)
    
    sum2 = np.sum(A*(Lij-Dl),axis=1)      
    El = np.array([[sum2[0]],[sum2[1]],[sum2[2]],[sum2[3]],[sum2[4]]])+np.dot(B,(Li0-Dl0))
    
    # calculate important parameters
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
    usv_betad_dot = np.array([u0*usv1_nu_dot[1],
                              u0*usv2_nu_dot[1],
                              u0*usv3_nu_dot[1],
                              u0*usv4_nu_dot[1],
                              u0*usv5_nu_dot[1]])/(usv_Ud**2)
    
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
        rc[j] = psicf_dot[j]-k2[j]*e_psi[j]
    
    # NN controller
    usv1_wd = np.array([Uc[0],rc[0]])
    usv2_wd = np.array([Uc[1],rc[1]])
    usv3_wd = np.array([Uc[2],rc[2]])
    usv4_wd = np.array([Uc[3],rc[3]])
    usv5_wd = np.array([Uc[4],rc[4]])
    usv1_Delta_T = usv1_Tc1-usv1_Tf
    usv2_Delta_T = usv2_Tc1-usv2_Tf
    usv3_Delta_T = usv3_Tc1-usv3_Tf
    usv4_Delta_T = usv4_Tc1-usv4_Tf
    usv5_Delta_T = usv5_Tc1-usv5_Tf
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
        usv1_What = np.zeros((9,2))
        usv2_What = np.zeros((9,2))
        usv3_What = np.zeros((9,2))
        usv4_What = np.zeros((9,2))
        usv5_What = np.zeros((9,2))
        
    nnc1 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    nnc2 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    nnc3 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    nnc4 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
    nnc5 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,4]),A_aux=np.diag([1,0.5]))
        
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)
    usv2_Lambda = nnc2.auxiliary_system(usv2_Lambda,usv2_Delta_T,ts)
    usv3_Lambda = nnc3.auxiliary_system(usv3_Lambda,usv3_Delta_T,ts)
    usv4_Lambda = nnc4.auxiliary_system(usv4_Lambda,usv4_Delta_T,ts)
    usv5_Lambda = nnc5.auxiliary_system(usv5_Lambda,usv5_Delta_T,ts)
    
    # input num = 3
    usv1_X = np.array([usv1_nu[0],usv1_nu[1],usv1_nu[2]])
    usv2_X = np.array([usv2_nu[0],usv2_nu[1],usv2_nu[2]])
    usv3_X = np.array([usv3_nu[0],usv3_nu[1],usv3_nu[2]])
    usv4_X = np.array([usv4_nu[0],usv4_nu[1],usv4_nu[2]])
    usv5_X = np.array([usv5_nu[0],usv5_nu[1],usv5_nu[2]]) # NN inputs
    '''
    # input num = 6  
    usv1_X = np.array([usv1_nu[0],usv1_nu[1],usv1_nu[2],usv1_nu[0]*usv1_nu[1],\
                       usv1_nu[0]*usv1_nu[2],usv1_nu[1]*usv1_nu[2]])
    usv2_X = np.array([usv2_nu[0],usv2_nu[1],usv2_nu[2],usv2_nu[0]*usv2_nu[1],\
                       usv2_nu[0]*usv2_nu[2],usv2_nu[1]*usv2_nu[2]])
    usv3_X = np.array([usv3_nu[0],usv3_nu[1],usv3_nu[2],usv3_nu[0]*usv3_nu[1],\
                       usv3_nu[0]*usv3_nu[2],usv3_nu[1]*usv3_nu[2]])
    usv4_X = np.array([usv4_nu[0],usv4_nu[1],usv4_nu[2],usv4_nu[0]*usv4_nu[1],\
                       usv4_nu[0]*usv4_nu[2],usv4_nu[1]*usv4_nu[2]])
    usv5_X = np.array([usv5_nu[0],usv5_nu[1],usv5_nu[2],usv5_nu[0]*usv5_nu[1],\
                       usv5_nu[0]*usv5_nu[2],usv5_nu[1]*usv5_nu[2]]) # NN inputs
    '''
    
    usv1_w = np.array([usv_U[0],usv1_nu[2]])
    usv2_w = np.array([usv_U[1],usv2_nu[2]])
    usv3_w = np.array([usv_U[2],usv3_nu[2]])
    usv4_w = np.array([usv_U[3],usv4_nu[2]])
    usv5_w = np.array([usv_U[4],usv5_nu[2]])  # w
    
    usv1_ew = usv1_w-usv1_wd-usv1_Lambda
    usv2_ew = usv2_w-usv2_wd-usv2_Lambda
    usv3_ew = usv3_w-usv3_wd-usv3_Lambda
    usv4_ew = usv4_w-usv4_wd-usv4_Lambda
    usv5_ew = usv5_w-usv5_wd-usv5_Lambda  # tracking errors
    
    usv1_wdf, usv1_wdf_dot = nnc1.TD(usv1_wd,usv1_wdf,ts)
    usv2_wdf, usv2_wdf_dot = nnc2.TD(usv2_wd,usv2_wdf,ts)
    usv3_wdf, usv3_wdf_dot = nnc3.TD(usv3_wd,usv3_wdf,ts)
    usv4_wdf, usv4_wdf_dot = nnc4.TD(usv4_wd,usv4_wdf,ts)
    usv5_wdf, usv5_wdf_dot = nnc5.TD(usv5_wd,usv5_wdf,ts)
    
    # estimated disturbances
    usv1_Fhat, usv1_What = nnc1.estimator(usv1_X,usv1_What,usv1_ew,ts)
    usv2_Fhat, usv2_What = nnc2.estimator(usv2_X,usv2_What,usv2_ew,ts)
    usv3_Fhat, usv3_What = nnc3.estimator(usv3_X,usv3_What,usv3_ew,ts)
    usv4_Fhat, usv4_What = nnc4.estimator(usv4_X,usv4_What,usv4_ew,ts)
    usv5_Fhat, usv5_What = nnc5.estimator(usv5_X,usv5_What,usv5_ew,ts)
    
    # real disturbances
    usv1_F = np.array([usv1_fU,usv1_fr])
    usv2_F = np.array([usv2_fU,usv2_fr])
    usv3_F = np.array([usv3_fU,usv3_fr])
    usv4_F = np.array([usv4_fU,usv4_fr])
    usv5_F = np.array([usv5_fU,usv5_fr])

    # controller
    usv1_Tc0 = nnc1.control_law(usv1_ew,usv1_F,usv1_wdf_dot,usv1_Lambda)
    usv2_Tc0 = nnc2.control_law(usv2_ew,usv2_F,usv2_wdf_dot,usv2_Lambda)
    usv3_Tc0 = nnc3.control_law(usv3_ew,usv3_F,usv3_wdf_dot,usv3_Lambda)
    usv4_Tc0 = nnc4.control_law(usv4_ew,usv4_F,usv4_wdf_dot,usv4_Lambda)
    usv5_Tc0 = nnc5.control_law(usv5_ew,usv5_F,usv5_wdf_dot,usv5_Lambda)
    
    # Event triggered mechanism
    Tc0 = np.array([[float(usv1_Tc0[0]),float(usv2_Tc0[0]),float(usv3_Tc0[0]),\
                        float(usv4_Tc0[0]),float(usv5_Tc0[0])],\
                       [float(usv1_Tc0[1]),float(usv2_Tc0[1]),float(usv3_Tc0[1]),\
                        float(usv4_Tc0[1]),float(usv5_Tc0[1])]])
        
    tauc0 = np.dot(B1_usv,Tc0)
    
    usv1_tauc0 = tauc0[:,0]
    usv2_tauc0 = tauc0[:,1]
    usv3_tauc0 = tauc0[:,2]
    usv4_tauc0 = tauc0[:,3]
    usv5_tauc0 = tauc0[:,4]
    
    if t == 0:
        tauc1 = np.zeros((2,5))
        triggeredtime1 = np.zeros((5,int(Ns)))
        triggeredtime2 = np.zeros((5,int(Ns)))
        d_coef = np.array([[0.5],[0.5]])
        d     = d_coef*np.ones((2,5))
        gamma_coef = np.array([[0.01],[0.01]])
        gamma = gamma_coef*np.ones((2,5))
        
    tauc0_1 = (1+gamma)*tauc0
           
    tauc1,triggeredtime1,triggeredtime2 = ETM1(tauc0_1,tauc1,triggeredtime1,triggeredtime2,d,gamma,time=i)
    
    Tc1 = np.dot(np.linalg.inv(B1_usv),tauc1)
    
    usv1_Tc1 = np.array([[Tc1[0,0]],[Tc1[1,0]]])
    usv2_Tc1 = np.array([[Tc1[0,1]],[Tc1[1,1]]])
    usv3_Tc1 = np.array([[Tc1[0,2]],[Tc1[1,2]]])
    usv4_Tc1 = np.array([[Tc1[0,3]],[Tc1[1,3]]])
    usv5_Tc1 = np.array([[Tc1[0,4]],[Tc1[1,4]]])
    
    usv1_tauc1 = tauc1[:,0]
    usv2_tauc1 = tauc1[:,1]
    usv3_tauc1 = tauc1[:,2]
    usv4_tauc1 = tauc1[:,3]
    usv5_tauc1 = tauc1[:,4]
    
    # store time series
    xout[i,:] = np.hstack((np.array([t]),x0[0],y0[0],u0[0]))
    usv1_out[i,:] = np.hstack((usv1_nu[0],usv1_nu[1],usv1_nu[2],\
                               usv1_eta[0],usv1_eta[1],usv1_eta[2],\
                               usv1_Tc0[0],usv1_Tc0[1],usv1_Tc1[0],\
                               usv1_Tc1[1],usv1_Tf[0],usv1_Tf[1],\
                               usv1_fU[0],usv1_fr[0],usv1_Fhat[0],\
                               usv1_Fhat[1]))
    usv2_out[i,:] = np.hstack((usv2_nu[0],usv2_nu[1],usv2_nu[2],\
                               usv2_eta[0],usv2_eta[1],usv2_eta[2],\
                               usv2_Tc0[0],usv2_Tc0[1],usv2_Tc1[0],\
                               usv2_Tc1[1],usv2_Tf[0],usv2_Tf[1],\
                               usv2_fU[0],usv2_fr[0],usv2_Fhat[0],\
                               usv2_Fhat[1]))
    usv3_out[i,:] = np.hstack((usv3_nu[0],usv3_nu[1],usv3_nu[2],\
                               usv3_eta[0],usv3_eta[1],usv3_eta[2],\
                               usv3_Tc0[0],usv3_Tc0[1],usv3_Tc1[0],\
                               usv3_Tc1[1],usv3_Tf[0],usv3_Tf[1],\
                               usv3_fU[0],usv3_fr[0],usv3_Fhat[0],\
                               usv3_Fhat[1]))
    usv4_out[i,:] = np.hstack((usv4_nu[0],usv4_nu[1],usv4_nu[2],\
                               usv4_eta[0],usv4_eta[1],usv4_eta[2],\
                               usv4_Tc0[0],usv4_Tc0[1],usv4_Tc1[0],\
                               usv4_Tc1[1],usv4_Tf[0],usv4_Tf[1],\
                               usv4_fU[0],usv4_fr[0],usv4_Fhat[0],\
                               usv4_Fhat[1]))
    usv5_out[i,:] = np.hstack((usv5_nu[0],usv5_nu[1],usv5_nu[2],\
                               usv5_eta[0],usv5_eta[1],usv5_eta[2],\
                               usv5_Tc0[0],usv5_Tc0[1],usv5_Tc1[0],\
                               usv5_Tc1[1],usv5_Tf[0],usv5_Tf[1],\
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
fig = plt.figure(figsize=(6,6),dpi = 600)
plt.plot(usv1_y,usv1_x,'r-',linewidth = linewid)
plt.plot(usv2_y,usv2_x,'g-',linewidth = linewid)
plt.plot(usv3_y,usv3_x,'b-',linewidth = linewid)
plt.plot(usv4_y,usv4_x,'c-',linewidth = linewid)
plt.plot(usv5_y,usv5_x,'m-',linewidth = linewid)
plt.plot(y0,x0,'k--',linewidth = linewid)

for j in np.arange(0,int(Ns),3000):
    plt.plot(usv1_y[j],usv1_x[j],'rs')
    plt.plot(usv2_y[j],usv2_x[j],'go')
    plt.plot(usv3_y[j],usv3_x[j],'b>')
    plt.plot(usv4_y[j],usv4_x[j],'cp')
    plt.plot(usv5_y[j],usv5_x[j],'md')
    
plt.plot(usv1_y[0],usv1_x[0],'rs',label="usv1")
plt.plot(usv2_y[0],usv2_x[0],'go',label="usv2")
plt.plot(usv3_y[0],usv3_x[0],'b>',label="usv3")
plt.plot(usv4_y[0],usv4_x[0],'cp',label="usv4")
plt.plot(usv5_y[0],usv5_x[0],'md',label="usv5")

plt.legend(loc="upper right")
plt.xlabel("east (m)",font)
plt.ylabel("north (m)",font)

#################### thrust
fig = plt.figure(figsize=(6,6),dpi = 600)
f, ax = plt.subplots(2,1,sharex = True)

plt.subplot(211)
ax[0].plot(t,usv1_Tc0_1,'r-',linewidth = linewid,label="Tc0_1")
ax[0].plot(t,usv1_Tc1_1,'g-',linewidth = linewid,label="Tc1_1")
plt.legend(loc="upper right")
plt.ylim(-80,180)
plt.xlim(0,tfinal)

# 绘制缩放图
axins = ax[0].inset_axes((0.4, 0.2, 0.4, 0.3))
axins.plot(t,usv1_Tc0_1,'r-',linewidth = linewid,label="Tc0_1")
axins.plot(t,usv1_Tc1_1,'g-',linewidth = linewid,label="Tc1_1")
axins.set_xlim(20, 30)
axins.set_ylim(85, 95)

plt.subplot(212)
ax[1].plot(t,usv1_Tc0_2,'b-',linewidth = linewid,label="Tc0_2")
ax[1].plot(t,usv1_Tc1_2,'m-',linewidth = linewid,label="Tc1_2")
plt.legend(loc="upper right")
plt.ylim(-80,180)
plt.xlim(0,tfinal)

# 绘制缩放图
axins = ax[1].inset_axes((0.4, 0.2, 0.4, 0.3))
axins.plot(t,usv1_Tc0_2,'r-',linewidth = linewid,label="Tc0_1")
axins.plot(t,usv1_Tc1_2,'g-',linewidth = linewid,label="Tc1_1")
axins.set_xlim(25, 30)
axins.set_ylim(90, 110)

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

# velocies
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,usv1_U,'r-',linewidth = linewid,label = "usv1")
axarr[0].plot(t,usv2_U,'g-',linewidth = linewid,label = "usv2")
axarr[0].plot(t,usv3_U,'b-',linewidth = linewid,label = "usv3")
axarr[0].plot(t,usv4_U,'c-',linewidth = linewid,label = "usv4")
axarr[0].plot(t,usv5_U,'m-',linewidth = linewid,label = "usv5")
#axarr[0].plot(t,u0,'k--',label = "U0",linewidth = linewid)
plt.legend(loc="upper right")
plt.ylabel("U (m/s)",font,labelpad=18)

plt.subplot(312)
axarr[1].plot(t,usv1_v,'r-',linewidth = linewid)
axarr[1].plot(t,usv2_v,'g-',linewidth = linewid)
axarr[1].plot(t,usv3_v,'b-',linewidth = linewid)
axarr[1].plot(t,usv4_v,'c-',linewidth = linewid)
axarr[1].plot(t,usv5_v,'m-',linewidth = linewid)

plt.ylabel("v (m/s)",font,labelpad=6)

plt.subplot(313)
axarr[2].plot(t,usv1_r,'r-',linewidth = linewid)
axarr[2].plot(t,usv2_r,'g-',linewidth = linewid)
axarr[2].plot(t,usv3_r,'b-',linewidth = linewid)
axarr[2].plot(t,usv4_r,'c-',linewidth = linewid)
axarr[2].plot(t,usv5_r,'m-',linewidth = linewid)

plt.ylabel("r (rad/s)",font,labelpad=18)

plt.xlabel("time (s)", font)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.8)

# coordinated following errors
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,ef1,'r-',linewidth = linewid,label="$e_{f1}$")
axarr[0].plot(t,ef2,'g-',linewidth = linewid,label="$e_{f2}$")
axarr[0].plot(t,ef3,'b-',linewidth = linewid,label="$e_{f3}$")
axarr[0].plot(t,ef4,'c-',linewidth = linewid,label="$e_{f4}$")
axarr[0].plot(t,ef5,'m-',linewidth = linewid,label="$e_{f5}$")

plt.ylabel("lati-error (m)",font,labelpad=18)
plt.legend(loc="upper right")

plt.subplot(312)
axarr[1].plot(t,el1,'r-',linewidth = linewid,label = "$e_{l1}$")
axarr[1].plot(t,el2,'g-',linewidth = linewid,label = "$e_{l2}$")
axarr[1].plot(t,el3,'b-',linewidth = linewid,label = "$e_{l3}$")
axarr[1].plot(t,el4,'c-',linewidth = linewid,label = "$e_{l4}$")
axarr[1].plot(t,el5,'m-',linewidth = linewid,label = "$e_{l5}$")

plt.ylabel("logi-error (m)",font,labelpad=18)
plt.legend(loc="upper right")

plt.subplot(313)
axarr[2].plot(t,chie1,'r-',linewidth = linewid,label = "$\chi_{e1}$")
axarr[2].plot(t,chie2,'g-',linewidth = linewid,label = "$\chi_{e2}$")
axarr[2].plot(t,chie3,'b-',linewidth = linewid,label = "$\chi_{e3}$")
axarr[2].plot(t,chie4,'c-',linewidth = linewid,label = "$\chi_{e4}$")
axarr[2].plot(t,chie5,'m-',linewidth = linewid,label = "$\chi_{e5}$")

plt.legend(loc="upper right")
plt.ylabel(r"$\theta_{e} (rad)$",font,labelpad=18)

plt.xlabel("time (s)",font,labelpad=18)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.8)

# estimation erros

fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(2,1,sharex = True)

plt.subplot(211)
axarr[0].plot(t,usv1_fU,'r-',linewidth = linewid,label = r"$f_{u}$")
axarr[0].plot(t,usv1_fUhat,'b--',linewidth = linewid,label=r"$\hat{f}_{u}$")
plt.ylabel("lati-error (m)",font,labelpad=18)
plt.title("usv1")
plt.legend(loc="upper right")
plt.xlim(0,tfinal)

plt.subplot(212)
axarr[1].plot(t,usv1_fr,'r-',linewidth = linewid,label=r"$f_{r}$")
axarr[1].plot(t,usv1_frhat,'b--',linewidth = linewid,label=r"$\hat{f}_{r}$")
plt.ylabel("lati-error (m)",font,labelpad=18)
plt.title("usv1")
plt.legend(loc="upper right")
plt.xlim(0,tfinal)

# triggered time 1
fig = plt.figure(figsize=(6,6),dpi = 600)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['r|','g|','b|','c|','m|'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime1[j1,j2]:
            plt.plot(t[j2],triggeredtime1[j1,j2],line_style[j1],\
                     markersize=5)
    plt.plot(t[j2],triggeredtime1[j1,j2],line_style[j1],\
             label = title[j1], markersize=8)
plt.legend(loc = "upper right")
plt.ylim(0,8)
plt.xlim(0,tfinal)

'''
# triggered time 2
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['rx','gx','bx','cx','mx'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    plt.subplot(5,1,(j1+1))
    for j2 in range(int(Ns)):
        if triggeredtime2[j1,j2]==1:
            axarr[j1].plot(Tk[j2],triggeredtime2[j1,j2],line_style[j1],markersize=2)
    plt.title(title[j1])
    plt.ylim(0,2)
    
f.subplots_adjust(hspace=0.8)
     
# triggered time 3   
triggeredtime3 = np.zeros((5,int(Ns)))
for j1 in range(5):
    for j2 in range(int(Ns)):
        if triggeredtime1[j1,j2]==1 or triggeredtime2[j1,j2]==1:
            triggeredtime3[j1,j2] = 1
            
fig = plt.figure(figsize=(6,6),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

title = np.array([r"$usv1$",r"$usv2$",r"$usv3$",\
                  r"$usv4$",r"$usv5$"])
line_style = np.array(['rx','gx','bx','cx','mx'])
Tk = np.arange(1,int(Ns+1),1)
for j1 in range(5):
    plt.subplot(5,1,(j1+1))
    for j2 in range(int(Ns)):
        if triggeredtime3[j1,j2]==1:
            axarr[j1].plot(Tk[j2],triggeredtime3[j1,j2],line_style[j1],markersize=2)
    plt.title(title[j1])
    plt.ylim(0,2)

'''

    
    
    
    
    
    
    
    
    
    
    
    
    
    