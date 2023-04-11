# -*- coding: utf-8 -*-
"""
Created on Sun Oct 30 11:38:52 2022

@author: qys
"""
import numpy as np


class ESOcontroller1:
    def __init__(self,Ko1,Ko2,KW,A_aux):
        # usv parameters
        self.m11 = 50.5
        self.dp = 0.26
        self.m33 = 17.21
        self.B = np.array([[1/self.m11,1/self.m11],\
                           [self.dp/self.m33,-self.dp/self.m33]])
        # ESO parameters
        self.Ko1 = Ko1
        self.Ko2 = Ko2
        # auxiliary system parameters
        self.A = A_aux
        # TD
        self.Twf = 0.1
        # controller parameters
        self.KW = KW
        
    def estimator(self,what,w,Fhat,T,ts):
        ew = w-what
        what_dot = Fhat+np.dot(self.Ko1,ew)+np.dot(self.B,T)
        what = ts*what_dot+what
        Fhat_dot = np.dot(self.Ko2,ew)
        Fhat = ts*Fhat_dot+Fhat
        return Fhat,what
   
    def auxiliary_system(self,Lambda,Delta_T,ts):
        Lambda_dot = -np.dot(self.A,Lambda)-np.dot(self.B,Delta_T)
        Lambda = ts*Lambda_dot+Lambda
        
        return Lambda
    
    def TD(self,wd,wdf,ts):
        wdf_dot = -(wdf-wd)/self.Twf
        wdf = ts*wdf_dot+wdf
        
        return wdf,wdf_dot
    
    def control_law(self,ew,Fhat,wd_dot,Lambda):
        Tc = np.dot(np.linalg.inv(self.B),(wd_dot-\
                    np.dot(self.A,Lambda)-Fhat-np.dot(self.KW,ew)))
        return Tc
        
    
    
    