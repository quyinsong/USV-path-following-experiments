# -*- coding: utf-8 -*-
"""
Created on Sun Oct 23 11:44:22 2022

@author: qys
"""
import numpy as np

class NNcontroller1:
     '''
     construct NNcontroller class
     '''
     
     def __init__(self,inputs_num,hiddenlayer_num,outputs_num,Gamma,kw,KW,A_aux):
         # usv parameters
         self.m11 = 50.5
         self.dp = 0.26
         self.m33 = 17.21
         self.B = np.array([[1/self.m11,1/self.m11],\
                            [self.dp/self.m33,-self.dp/self.m33]])
         # network parameters
         self.inputs_number = inputs_num
         self.hiddenlayer_number = hiddenlayer_num
         self.outputs_number = outputs_num
         self.Gamma = Gamma
         self.kw = kw
         self.c = np.zeros((self.inputs_number,self.hiddenlayer_number))
         self.b = 2*np.ones((self.hiddenlayer_number,1))
         # auxiliary system parameters
         self.A = A_aux
         # TD
         self.Twf = 0.1
         # controller parameters
         self.KW = KW
         
     def estimator(self,X,What,ew,ts):
        # estimate
        PhiX = np.zeros((self.hiddenlayer_number,1))
        for i in range(self.hiddenlayer_number):
            PhiX[i] = np.exp(-(np.power(np.linalg.norm((X-self.c[:,i]),ord=2),2))/np.power(self.b[i],2))
        Fhat = np.dot(What.T,PhiX)
        # weights update
        What_dot = self.Gamma*(np.dot(PhiX,ew.T)-self.kw*What)
        What = ts*What_dot+What
        # return estimate value and next step weights
        return Fhat,What
    
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
  
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         
         