# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 19:55:17 2022

@author: qys
"""
# 在CPF_test1 基础上增加观测器，估计侧滑引起的不确定项
# 结论： 效果不好，容易引起振荡，且没有必要估计
import numpy as np
import matplotlib.pyplot as plt
from usvmodel import USV1
from nncontroller import NNcontroller1

# initialize
ts = 0.02
tfinal = 100
Ns = tfinal/ts

Nxout = 18
Nusv1 = 14
Nusv2 = 14
Nusv3 = 14
Nusv4 = 14
Nusv5 = 14
Nerror = 15
xout = np.zeros((int(Ns),Nxout))
usv1_out = np.zeros((int(Ns),Nusv1))
usv2_out = np.zeros((int(Ns),Nusv2))
usv3_out = np.zeros((int(Ns),Nusv3))
usv4_out = np.zeros((int(Ns),Nusv4))
usv5_out = np.zeros((int(Ns),Nusv5))
error_out = np.zeros((int(Ns),Nerror))

# path information
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
        Sum += h*(fx3_integrated(xi+h*i)+fx3_integrated(xi+h*(i+1)))/2
    
    return Sum               # return arc length
    

# simulation start
for i in range(int(Ns)):
    t = i*ts
    # initializing
    if t == 0:
        usv1_nu = np.array([[0],[0],[0]]) # usv1-usv5 states initialization
        usv1_eta = np.array([[5],[0],[0]])
        usv1_Tf = np.array([[0],[0]])
        usv1_Tc = np.array([[0],[0]])
        usv2_nu = np.array([[0],[0],[0]])
        usv2_eta = np.array([[0],[-15],[0]])
        usv2_Tf = np.array([[0],[0]])
        usv2_Tc = np.array([[0],[0]])
        usv3_nu = np.array([[0],[0],[0]])
        usv3_eta = np.array([[0],[-30],[0]])
        usv3_Tf = np.array([[0],[0]])
        usv3_Tc = np.array([[0],[0]])
        usv4_nu = np.array([[0],[0],[0]])
        usv4_eta = np.array([[0],[20],[0]])
        usv4_Tf = np.array([[0],[0]])
        usv4_Tc = np.array([[0],[0]])
        usv5_nu = np.array([[0],[0],[0]])
        usv5_eta = np.array([[0],[30],[0]])
        usv5_Tf = np.array([[0],[0]])
        usv5_Tc = np.array([[0],[0]])
        x0 = np.array([0])  # virtual target position initialization
        y0 = np.array([10])
        u0 = np.array([0.5])
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

        Df = np.array([[0,0,0,0,0],
                       [-8,0,0,0,0],
                       [0,-8,0,0,0],
                       [8,0,0,0,0],
                       [0,0,0,8,0]])  # desired lati-error
        Df0 = np.array([[0],[0],[0],[0],[0]])
        Dl = np.array([[0,0,0,0,0],
                       [-8,0,0,0,0],
                       [0,-8,0,0,0],
                       [-8,0,0,0,0],
                       [0,0,0,-8,0]]) # desired logi-error
        Dl0 = np.array([[0],[0],[0],[0],[0]])
        
        k1 = np.array([[0.2],[0.2],[0.2],[0.2],[0.2]]) # guidane law parameters
        k2 = np.array([[0.5],[0.5],[0.5],[0.5],[0.5]])
        k3 = np.array([[1],[1],[1],[1],[1]])
        
    # disturbancs and currents
    
    tw = np.array([[15*np.sin(0.1*t)*np.cos(0.05*t)+10],
                   [10*np.sin(0.1*t)*np.cos(0.05*t)+10],
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
    usv1_nu,usv1_eta,usv1_nu_dot,usv1_eta_dot,usv1_Tf,usv1_fu,usv1_fr = \
        usv1.state_update1(usv1_nu,usv1_eta,usv1_Tf,usv1_Tc,vc,tw,ts)
    usv2_nu,usv2_eta,usv2_nu_dot,usv2_eta_dot,usv2_Tf,usv2_fu,usv2_fr = \
        usv2.state_update1(usv2_nu,usv2_eta,usv2_Tf,usv2_Tc,vc,tw,ts)
    usv3_nu,usv3_eta,usv3_nu_dot,usv3_eta_dot,usv3_Tf,usv3_fu,usv3_fr = \
        usv3.state_update1(usv3_nu,usv3_eta,usv3_Tf,usv3_Tc,vc,tw,ts)
    usv4_nu,usv4_eta,usv4_nu_dot,usv4_eta_dot,usv4_Tf,usv4_fu,usv4_fr = \
        usv4.state_update1(usv4_nu,usv4_eta,usv4_Tf,usv4_Tc,vc,tw,ts)
    usv5_nu,usv5_eta,usv5_nu_dot,usv5_eta_dot,usv5_Tf,usv5_fu,usv5_fr = \
        usv5.state_update1(usv5_nu,usv5_eta,usv5_Tf,usv5_Tc,vc,tw,ts)
    
    # virtual target state update
    f0,fx0,fy0,ffxx0,ffyy0,ffxy0 = fxy3(x0,y0)
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
        F[j],usv_fx[j],usv_fy[j],usv_ffxx[j],usv_ffyy[j],usv_ffxy[j] = fxy3(eta[0,j],eta[1,j])
    
    sum1 = np.sum(A*Df,axis=1)
    Ef = np.dot((L+B),F)-np.array([[sum1[0]],[sum1[1]],[sum1[2]],[sum1[3]],[sum1[4]]])-np.dot(B,Df0)
            
    # calculate longi-coordinated errors
    Lij = np.zeros((5,5))   # request memory to store logi-coordinated errors
    Li0 = np.zeros((5,1))
    
    for m in range(5):
        for n in range(5):
            Lij[m,n] = arc_length(eta[0,m],eta[0,n],20)
            Li0[m,0] = arc_length(eta[0,m],x0,20)
    
    sum2 = np.sum(A*(Lij-Dl),axis=1)      
    El = np.array([[sum2[0]],[sum2[1]],[sum2[2]],[sum2[3]],[sum2[4]]])+np.dot(B,(Li0-Dl0))
    
    # calculate important parameters
    chid = np.zeros((5,1))
    chid_dot = np.zeros((5,1))
    fxynorm2 = np.zeros((5,1))
    betad = np.zeros((5,1))
    chie = np.zeros((5,1))
    F_dot = np.zeros((5,1))
    
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
        betad[j] = np.arctan(nu[1,j]/u0)
        chie[j] = eta[2,j]-chid[j]
    
    # extended state observer
    if t == 0:
        Fhat = F
        thetahat = np.zeros((5,1))
        thetahat_dot = np.zeros((5,1))
        Fhat_dot = np.zeros((5,1))
        theta = np.zeros((5,1))
        ko1 = np.array([[2.]])
        ko2 = np.array([[20.]])
        
    for j in range(5):
        Fhat_dot[j] = np.sqrt(fxynorm2[j])*nu[0,j]*np.sin(chie[j])+\
            np.sqrt(fxynorm2[j])*thetahat[j]+ko1[0]*(F[j]-Fhat[j])
        thetahat_dot[j] = ko2[0]*(F[j]-Fhat[j])
        Fhat[j] = ts*Fhat_dot[j]+Fhat[j]
        thetahat[j] = ts*thetahat_dot[j]+thetahat[j]
        theta[j] = nu[1,j]*np.cos(chie[j])
        if thetahat[j]>=u0:
            thetahat[j] = 0
        
    # coordinated guidance law
    if t == 0:
        sigma = np.zeros((5,1))
        sigmaf = np.zeros((5,1))
        sigmaf_dot = np.zeros((5,1))
        Tf = np.array([0.1])
        
    for j in range(5):
        F_dot[j] = np.sqrt(fxynorm2[j])*(nu[0,j]*np.sin(chie[j])+thetahat[j])
        sigma[j] = (k1[j]*F[j]*(thetahat[j]**2)+np.sqrt((u0**2)-(thetahat[j]**2)
                    +(k1[j]*F[j]*u0)**2)*np.abs(thetahat[j]))/(u0**2-thetahat[j]**2)
        sigmaf_dot[j] = -(sigmaf[j]-sigma[j])/Tf
        sigmaf[j] = ts*sigmaf_dot[j]+sigmaf[j]
        
    Ef_dot = np.dot((L+B),F_dot)
    
    uc = np.zeros((5,1))
    rc = np.zeros((5,1))
    s = np.zeros((5,1))
    
    for j in range(5):
        s[j] = np.arctan(k1[j]*Ef[j]+sigma[j])+chie[j]
        uc[j] = u0*(1-np.tanh(k3[j]*El[j]))
        rc[j] = -(k1[j]*Ef_dot[j])/(1+np.square(k1[j]*Ef[j]+sigma[j]))-\
                k2[j]*np.sign(s[j])*(np.abs(s[j])**1)+chid_dot[j]
    
    # NN controller
    usv1_wd = np.array([uc[0],rc[0]])
    usv2_wd = np.array([uc[1],rc[1]])
    usv3_wd = np.array([uc[2],rc[2]])
    usv4_wd = np.array([uc[3],rc[3]])
    usv5_wd = np.array([uc[4],rc[4]])
    usv1_Delta_T = usv1_Tc-usv1_Tf
    usv2_Delta_T = usv2_Tc-usv2_Tf
    usv3_Delta_T = usv3_Tc-usv3_Tf
    usv4_Delta_T = usv4_Tc-usv4_Tf
    usv5_Delta_T = usv5_Tc-usv5_Tf
    if t==0:
        B_usv = np.array([[1/50.5,1/50.5],[0.26/17.21,-0.26/17.21]])
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
                            KW=np.diag([4,8]),A_aux=np.diag([1,0.5]))
    nnc2 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,8]),A_aux=np.diag([1,0.5]))
    nnc3 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,8]),A_aux=np.diag([1,0.5]))
    nnc4 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,8]),A_aux=np.diag([1,0.5]))
    nnc5 = NNcontroller1(inputs_num=3,hiddenlayer_num=9,\
                            outputs_num=2,Gamma=2.5,kw=0.005,\
                            KW=np.diag([4,8]),A_aux=np.diag([1,0.5]))
        
    usv1_Lambda = nnc1.auxiliary_system(usv1_Lambda,usv1_Delta_T,ts)
    usv2_Lambda = nnc2.auxiliary_system(usv2_Lambda,usv2_Delta_T,ts)
    usv3_Lambda = nnc3.auxiliary_system(usv3_Lambda,usv3_Delta_T,ts)
    usv4_Lambda = nnc4.auxiliary_system(usv4_Lambda,usv4_Delta_T,ts)
    usv5_Lambda = nnc5.auxiliary_system(usv5_Lambda,usv5_Delta_T,ts)
    
    usv1_X = np.array([usv1_nu[0],usv1_nu[2],usv1_nu[0]*usv1_nu[2]])
    usv2_X = np.array([usv2_nu[0],usv2_nu[2],usv2_nu[0]*usv2_nu[2]])
    usv3_X = np.array([usv3_nu[0],usv3_nu[2],usv3_nu[0]*usv3_nu[2]])
    usv4_X = np.array([usv4_nu[0],usv4_nu[2],usv4_nu[0]*usv4_nu[2]])
    usv5_X = np.array([usv5_nu[0],usv5_nu[2],usv5_nu[0]*usv5_nu[2]]) # NN inputs
    
    usv1_w = np.array([usv1_nu[0],usv1_nu[2]])
    usv2_w = np.array([usv2_nu[0],usv2_nu[2]])
    usv3_w = np.array([usv3_nu[0],usv3_nu[2]])
    usv4_w = np.array([usv4_nu[0],usv4_nu[2]])
    usv5_w = np.array([usv5_nu[0],usv5_nu[2]])  # w
    
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
    
    usv1_Fhat, usv1_What = nnc1.estimator(usv1_X,usv1_What,usv1_ew,ts)
    usv2_Fhat, usv2_What = nnc2.estimator(usv2_X,usv2_What,usv2_ew,ts)
    usv3_Fhat, usv3_What = nnc3.estimator(usv3_X,usv3_What,usv3_ew,ts)
    usv4_Fhat, usv4_What = nnc4.estimator(usv4_X,usv4_What,usv4_ew,ts)
    usv5_Fhat, usv5_What = nnc5.estimator(usv5_X,usv5_What,usv5_ew,ts)
    
    usv1_Tc = nnc1.control_law(usv1_ew,usv1_Fhat,usv1_wdf_dot,usv1_Lambda)
    usv2_Tc = nnc2.control_law(usv2_ew,usv2_Fhat,usv2_wdf_dot,usv2_Lambda)
    usv3_Tc = nnc3.control_law(usv3_ew,usv3_Fhat,usv3_wdf_dot,usv3_Lambda)
    usv4_Tc = nnc4.control_law(usv4_ew,usv4_Fhat,usv4_wdf_dot,usv4_Lambda)
    usv5_Tc = nnc5.control_law(usv5_ew,usv5_Fhat,usv5_wdf_dot,usv5_Lambda)
    
    # store time series
    xout[i,:] = np.hstack((np.array([t]),x0[0],y0[0],theta[0],theta[1],\
                           theta[2],theta[3],theta[4],thetahat[0],thetahat[1],\
                           thetahat[2],thetahat[3],thetahat[4],sigmaf_dot[0],\
                           sigmaf_dot[1],sigmaf_dot[2],sigmaf_dot[3],sigmaf_dot[4]))
    usv1_out[i,:] = np.hstack((usv1_nu[0],usv1_nu[1],usv1_nu[2],\
                               usv1_eta[0],usv1_eta[1],usv1_eta[2],\
                               usv1_Tc[0],usv1_Tc[1],usv1_Tf[0],\
                               usv1_Tf[1],usv1_fu[0],usv1_fr[0],\
                               usv1_Fhat[0],usv1_Fhat[1]))
    usv2_out[i,:] = np.hstack((usv2_nu[0],usv2_nu[1],usv2_nu[2],\
                               usv2_eta[0],usv2_eta[1],usv2_eta[2],\
                               usv2_Tc[0],usv2_Tc[1],usv2_Tf[0],\
                               usv2_Tf[1],usv2_fu[0],usv2_fr[0],\
                               usv2_Fhat[0],usv2_Fhat[1]))
    usv3_out[i,:] = np.hstack((usv3_nu[0],usv3_nu[1],usv3_nu[2],\
                               usv3_eta[0],usv3_eta[1],usv3_eta[2],\
                               usv3_Tc[0],usv3_Tc[1],usv3_Tf[0],\
                               usv3_Tf[1],usv3_fu[0],usv3_fr[0],\
                               usv3_Fhat[0],usv3_Fhat[1]))
    usv4_out[i,:] = np.hstack((usv4_nu[0],usv4_nu[1],usv4_nu[2],\
                               usv4_eta[0],usv4_eta[1],usv4_eta[2],\
                               usv4_Tc[0],usv4_Tc[1],usv4_Tf[0],\
                               usv4_Tf[1],usv4_fu[0],usv4_fr[0],\
                               usv4_Fhat[0],usv4_Fhat[1]))
    usv5_out[i,:] = np.hstack((usv5_nu[0],usv5_nu[1],usv5_nu[2],\
                               usv5_eta[0],usv5_eta[1],usv5_eta[2],\
                               usv5_Tc[0],usv5_Tc[1],usv5_Tf[0],\
                               usv5_Tf[1],usv5_fu[0],usv5_fr[0],\
                               usv5_Fhat[0],usv5_Fhat[1]))
    error_out[i,:] = np.hstack((Ef[0],Ef[1],Ef[2],Ef[3],Ef[4],\
                                El[0],El[1],El[2],El[3],El[4],\
                                s[0],s[1],s[2],s[3],s[4]))
    
    
    
t = xout[:,0]
x0 = xout[:,1]
y0 = xout[:,2]
theta1 = xout[:,3]
theta2 = xout[:,4]
theta3 = xout[:,5]
theta4 = xout[:,6]
theta5 = xout[:,7]
thetahat1 = xout[:,8]
thetahat2 = xout[:,9]
thetahat3 = xout[:,10]
thetahat4 = xout[:,11]
thetahat5 = xout[:,12]
sigmaf_dot1 = xout[:,13]
sigmaf_dot2 = xout[:,14]
sigmaf_dot3 = xout[:,15]
sigmaf_dot4 = xout[:,16]
sigmaf_dot5 = xout[:,17]

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

usv2_u = usv2_out[:,0]
usv2_v = usv2_out[:,1]
usv2_r = usv2_out[:,2]
usv2_x = usv2_out[:,3]
usv2_y = usv2_out[:,4]
usv2_psi = usv2_out[:,5]
usv2_Tc1 = usv2_out[:,6]
usv2_Tc2 = usv2_out[:,7]
usv2_Tf1 = usv2_out[:,8]
usv2_Tf2 = usv2_out[:,9]
usv2_fu = usv2_out[:,10]
usv2_fr = usv2_out[:,11]
usv2_fuhat = usv2_out[:,12]
usv2_frhat = usv2_out[:,13]

usv3_u = usv3_out[:,0]
usv3_v = usv3_out[:,1]
usv3_r = usv3_out[:,2]
usv3_x = usv3_out[:,3]
usv3_y = usv3_out[:,4]
usv3_psi = usv3_out[:,5]
usv3_Tc1 = usv3_out[:,6]
usv3_Tc2 = usv3_out[:,7]
usv3_Tf1 = usv3_out[:,8]
usv3_Tf2 = usv3_out[:,9]
usv3_fu = usv3_out[:,10]
usv3_fr = usv3_out[:,11]
usv3_fuhat = usv3_out[:,12]
usv3_frhat = usv3_out[:,13]

usv4_u = usv4_out[:,0]
usv4_v = usv4_out[:,1]
usv4_r = usv4_out[:,2]
usv4_x = usv4_out[:,3]
usv4_y = usv4_out[:,4]
usv4_psi = usv4_out[:,5]
usv4_Tc1 = usv4_out[:,6]
usv4_Tc2 = usv4_out[:,7]
usv4_Tf1 = usv4_out[:,8]
usv4_Tf2 = usv4_out[:,9]
usv4_fu = usv4_out[:,10]
usv4_fr = usv4_out[:,11]
usv4_fuhat = usv4_out[:,12]
usv4_frhat = usv4_out[:,13]

usv5_u = usv5_out[:,0]
usv5_v = usv5_out[:,1]
usv5_r = usv5_out[:,2]
usv5_x = usv5_out[:,3]
usv5_y = usv5_out[:,4]
usv5_psi = usv5_out[:,5]
usv5_Tc1 = usv5_out[:,6]
usv5_Tc2 = usv5_out[:,7]
usv5_Tf1 = usv5_out[:,8]
usv5_Tf2 = usv5_out[:,9]
usv5_fu = usv5_out[:,10]
usv5_fr = usv5_out[:,11]
usv5_fuhat = usv5_out[:,12]
usv5_frhat = usv5_out[:,13]

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
s1 = error_out[:,10]
s2 = error_out[:,11]
s3 = error_out[:,12]
s4 = error_out[:,13]
s5 = error_out[:,14]

# plots
font = {'family' : 'Times New Roman',
'weight' : 'normal',
'size'   : 12,
}

# north-east map
fig = plt.figure(figsize=(6,6
                          ),dpi = 600)
plt.plot(usv1_y,usv1_x,'r-')
plt.plot(usv2_y,usv2_x,'g-')
plt.plot(usv3_y,usv3_x,'b-')
plt.plot(usv4_y,usv4_x,'c-')
plt.plot(usv5_y,usv5_x,'m-')
plt.plot(y0,x0,'k--')

for j in np.arange(0,int(Ns),8000):
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

# thrust
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

plt.subplot(511)
axarr[0].plot(t,usv1_Tc1,label="Tc1")
axarr[0].plot(t,usv1_Tc2,label="Tc2")
axarr[0].plot(t,usv1_Tf1,label="Tf1")
axarr[0].plot(t,usv1_Tf2,label="Tf2")
plt.title("usv1")
plt.ylabel("thrust (N)")
plt.legend(loc="upper right")

plt.subplot(512)
axarr[1].plot(t,usv2_Tc1)
axarr[1].plot(t,usv2_Tc2)
axarr[1].plot(t,usv2_Tf1)
axarr[1].plot(t,usv2_Tf2)
plt.title("usv2")
plt.ylabel("thrust (N)")

plt.subplot(513)
axarr[2].plot(t,usv3_Tc1)
axarr[2].plot(t,usv3_Tc2)
axarr[2].plot(t,usv3_Tf1)
axarr[2].plot(t,usv3_Tf2)
plt.title("usv3")

plt.subplot(514)
axarr[3].plot(t,usv4_Tc1)
axarr[3].plot(t,usv4_Tc2)
axarr[3].plot(t,usv4_Tf1)
axarr[3].plot(t,usv4_Tf2)
plt.title("usv4")

plt.subplot(515)
axarr[4].plot(t,usv5_Tc1)
axarr[4].plot(t,usv5_Tc2)
axarr[4].plot(t,usv5_Tf1)
axarr[4].plot(t,usv5_Tf2)

plt.title("usv5")
plt.ylabel("thrust (N)")
plt.xlabel("time (s)")
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.8)

# velocies
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,usv1_u,'r-',label = "usv1")
axarr[0].plot(t,usv2_u,'g-',label = "usv2")
axarr[0].plot(t,usv3_u,'b-',label = "usv3")
axarr[0].plot(t,usv4_u,'c-',label = "usv4")
axarr[0].plot(t,usv5_u,'m-',label = "usv5")

plt.ylabel("u (m/s)",font,labelpad=18)

plt.subplot(312)
axarr[1].plot(t,usv1_v,'r-')
axarr[1].plot(t,usv2_v,'g-')
axarr[1].plot(t,usv3_v,'b-')
axarr[1].plot(t,usv4_v,'c-')
axarr[1].plot(t,usv5_v,'m-')

plt.ylabel("v (m/s)",font,labelpad=6)

plt.subplot(313)
axarr[2].plot(t,usv1_r,'r-')
axarr[2].plot(t,usv2_r,'g-')
axarr[2].plot(t,usv3_r,'b-')
axarr[2].plot(t,usv4_r,'c-')
axarr[2].plot(t,usv5_r,'m-')

plt.ylabel("r (rad/s)",font,labelpad=18)

plt.xlabel("time (s)", font)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.8)

# coordinated following errors
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(3,1,sharex = True)

plt.subplot(311)
axarr[0].plot(t,ef1,'r-',label="$e_{f1}$")
axarr[0].plot(t,ef2,'g-',label="$e_{f2}$")
axarr[0].plot(t,ef3,'b-',label="$e_{f3}$")
axarr[0].plot(t,ef4,'c-',label="$e_{f4}$")
axarr[0].plot(t,ef5,'m-',label="$e_{f5}$")

plt.ylabel("lati-error (m)",font,labelpad=18)
plt.legend(loc="upper right")

plt.subplot(312)
axarr[1].plot(t,el1,'r-',label = "$e_{l1}$")
axarr[1].plot(t,el2,'g-',label = "$e_{l2}$")
axarr[1].plot(t,el3,'b-',label = "$e_{l3}$")
axarr[1].plot(t,el4,'c-',label = "$e_{l4}$")
axarr[1].plot(t,el5,'m-',label = "$e_{l5}$")

plt.ylabel("logi-error (m)",font,labelpad=18)
plt.legend(loc="upper right")

plt.subplot(313)
axarr[2].plot(t,s1,'r-',label = "$s_{1}$")
axarr[2].plot(t,s2,'g-',label = "$s_{2}$")
axarr[2].plot(t,s3,'b-',label = "$s_{3}$")
axarr[2].plot(t,s4,'c-',label = "$s_{4}$")
axarr[2].plot(t,s5,'m-',label = "$s_{5}$")

plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.xlabel("time (s)",font,labelpad=18)
plt.xlim(0,tfinal)
f.subplots_adjust(hspace=0.8)

# NN estimation

fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

plt.subplot(511)
axarr[0].plot(t,usv1_fu,'r-',label = "$f_{u}$")
axarr[0].plot(t,usv1_fuhat,'b--',label="$\hat{f}_{u}$")
axarr[0].plot(t,usv1_fr,'g-',label="$f_{r}$")
axarr[0].plot(t,usv1_frhat,'m--',label="$\hat{f}_{r}$")
plt.ylabel("lati-error (m)",font,labelpad=18)
plt.title("usv1")
plt.legend(loc="upper right")

plt.subplot(512)
axarr[1].plot(t,usv2_fu,'r-',label = "$f_{u}$")
axarr[1].plot(t,usv2_fuhat,'b--',label="$\hat{f}_{u}$")
axarr[1].plot(t,usv2_fr,'g-',label="$f_{r}$")
axarr[1].plot(t,usv2_frhat,'m--',label="$\hat{f}_{r}$")
plt.title("usv2")
plt.legend(loc="upper right")

plt.subplot(513)
axarr[2].plot(t,usv3_fu,'r-',label = "$f_{u}$")
axarr[2].plot(t,usv3_fuhat,'b--',label="$\hat{f}_{u}$")
axarr[2].plot(t,usv3_fr,'g-',label="$f_{r}$")
axarr[2].plot(t,usv3_frhat,'m--',label="$\hat{f}_{r}$")
plt.title("usv3")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.subplot(514)
axarr[3].plot(t,usv4_fu,'r-',label = "$f_{u}$")
axarr[3].plot(t,usv4_fuhat,'b--',label="$\hat{f}_{u}$")
axarr[3].plot(t,usv4_fr,'g-',label="$f_{r}$")
axarr[3].plot(t,usv4_frhat,'m--',label="$\hat{f}_{r}$")
plt.title("usv4")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.subplot(515)
axarr[4].plot(t,usv5_fu,'r-',label = "$f_{u}$")
axarr[4].plot(t,usv5_fuhat,'b--',label="$\hat{f}_{u}$")
axarr[4].plot(t,usv5_fr,'g-',label="$f_{r}$")
axarr[4].plot(t,usv5_frhat,'m--',label="$\hat{f}_{r}$")
plt.title("usv5")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.xlabel("time (s)",font,labelpad=18)
plt.xlim(0,tfinal)

f.subplots_adjust(hspace=0.8)

# extended state observer
fig = plt.figure(figsize=(6,3),dpi = 600)
f, axarr = plt.subplots(5,1,sharex = True)

plt.subplot(511)
axarr[0].plot(t,theta1,'r-',label = r"$\theta_{1}$")
axarr[0].plot(t,thetahat1,'b--',label= r"$\hat{\theta}_{1}$")
plt.ylabel(r"theta (m/s)",font,labelpad=18)
plt.title("usv1")
plt.legend(loc="upper right")

plt.subplot(512)
axarr[1].plot(t,theta2,'r-',label = r"$\theta_{2}$")
axarr[1].plot(t,thetahat2,'b--',label= r"$\hat{\theta}_{2}$")

plt.title("usv2")
plt.legend(loc="upper right")

plt.subplot(513)
axarr[2].plot(t,theta3,'r-',label = r"$\theta_{3}$")
axarr[2].plot(t,thetahat3,'b--',label= r"$\hat{\theta}_{3}$")
plt.title("usv3")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.subplot(514)
axarr[3].plot(t,theta4,'r-',label = r"$\theta_{4}$")
axarr[3].plot(t,thetahat4,'b--',label= r"$\hat{\theta}_{4}$")
plt.title("usv4")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.subplot(515)
axarr[4].plot(t,theta1,'r-',label = r"$\theta_{5}$")
axarr[4].plot(t,thetahat1,'b--',label= r"$\hat{\theta}_{5}$")
plt.title("usv5")
plt.legend(loc="upper right")
plt.ylabel("\theta_{e} (rad)",font,labelpad=18)

plt.xlabel("time (s)",font,labelpad=18)
plt.xlim(0,tfinal)

f.subplots_adjust(hspace=0.8)
    
    
fig = plt.figure(figsize=(6,3),dpi = 600)
plt.plot(t,sigmaf_dot1,'r-',label = r"\dot{sigma}_{f1}") 
plt.plot(t,sigmaf_dot2,'g-',label = r"\dot{sigma}_{f1}") 
plt.plot(t,sigmaf_dot3,'b-',label = r"\dot{sigma}_{f1}") 
plt.plot(t,sigmaf_dot4,'c-',label = r"\dot{sigma}_{f1}") 
plt.plot(t,sigmaf_dot5,'m-',label = r"\dot{sigma}_{f1}") 

    
    
    
    
    
    
    
    
    
    
