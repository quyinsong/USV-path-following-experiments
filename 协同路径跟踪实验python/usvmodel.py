# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 14:10:28 2022

@author: qys
"""

import numpy as np

class USV1:
    '''
    construct the class of USV
    '''
    
    def __init__(self):
        '''
        

        Parameters
        ----------
        nu0 : vector[float]
            velocity vector: nu0 = [u0,v0,r0].T
        eta0 : vector[float]
            position vector: eta0 = [x0,y0,psi0].T

        Returns
        -------
        None.

        '''
        self.M = np.diag([50.5,84.36,17.21])      # inertial matrix
        self.D = np.diag([151.57,132.5,34.56])    # damping matrix
        self.u_max = 1.5                          # the maximum surge speed
        self.dp = 0.26
        self.B = np.array([[1,1],[0,0],[self.dp,-self.dp]])
        self.T_max = self.u_max*self.D[0,0]       # maximum thrust
        self.T_min = -55                          # minimum thrust
        self.Tf_time = 0.05                       # integration constant of thruster
        self.Tdot_max = 50                        # maximum change rate of thrust
        self.Tdot_min = -30                       # minimum change rate of thrust
        
    def state_update(self,nu,eta,Tf,Tc,vc,tw,ts):
        '''
        

        Parameters
        ----------
        nu : TYPE
            velocity state of USV
        eta : TYPE
            position state of USV
        Tf : TYPE
            thrust filtered value
        Tc : TYPE
            thrust command input
        vc : TYPE
            currents
        tw : TYPE
            unknown disturbances
        ts : TYPE
            simulation step

        Returns
        -------
        next_nu
        next_eta
        next_Tf
        fu
        fr

        '''
        Tf_dot = -(Tf-Tc)/self.Tf_time
        
        # limit the change rate of thrust command
        for i in range(2):
            if Tf_dot[i] >= self.Tdot_max:
                Tf_dot[i] = self.Tdot_max
            if Tf_dot[i] <= self.Tdot_min:
                Tf_dot[i] = self.Tdot_min
                
        # update the thrust
        Tf = ts*Tf_dot+Tf
        
        # limit the amplitude of thrust command
        for i in range(2):
            if Tf[i] >= self.T_max:
                Tf[i] = self.T_max
            if Tf[i] <= self.T_min:
                Tf[i] = self.T_min
                
        # USV's state update
        psi = eta[2]
        u = nu[0]
        v = nu[1]
        R = np.array([[np.cos(psi),-np.sin(psi),0],[np.sin(psi),np.cos(psi),0],[0,0,1]],dtype = 'float')        # rotation matrix
        C = np.array([[0,0,self.M[2,2]*v],[0,0,self.M[1,1]*u],[self.M[2,2]*v,-self.M[1,1]*u,0]],dtype = 'float')
        F = np.dot(np.linalg.inv(self.M),(-np.dot(self.D,nu)-np.dot(C,nu)+tw))
        fu = F[0]
        fr = F[2]
        
        nu_dot = F+np.dot(np.linalg.inv(self.M),np.dot(self.B,Tf))
        eta_dot = np.dot(R,nu)+vc
        
        nu = ts*nu_dot+nu
        eta = ts*eta_dot+eta
        
        return nu,eta,Tf,fu,fr
    
    def state_update1(self,nu,eta,Tf,Tc,vc,tw,ts):
        '''
        

        Parameters
        ----------
        nu : TYPE
            velocity state of USV
        eta : TYPE
            position state of USV
        Tf : TYPE
            thrust filtered value
        Tc : TYPE
            thrust command input
        vc : TYPE
            currents
        tw : TYPE
            unknown disturbances
        ts : TYPE
            simulation step

        Returns
        -------
        next_nu
        next_eta
        nu_dot
        eta_dot
        next_Tf
        fu
        fr

        '''
        Tf_dot = -(Tf-Tc)/self.Tf_time
        
        # limit the change rate of thrust command
        for i in range(2):
            if Tf_dot[i] >= self.Tdot_max:
                Tf_dot[i] = self.Tdot_max
            if Tf_dot[i] <= self.Tdot_min:
                Tf_dot[i] = self.Tdot_min
                
        # update the thrust
        Tf = ts*Tf_dot+Tf
        
        # limit the amplitude of thrust command
        for i in range(2):
            if Tf[i] >= self.T_max:
                Tf[i] = self.T_max
            if Tf[i] <= self.T_min:
                Tf[i] = self.T_min
                
        # USV's state update
        psi = eta[2]
        u = nu[0]
        v = nu[1]
        R = np.array([[np.cos(psi),-np.sin(psi),0],[np.sin(psi),np.cos(psi),0],[0,0,1]],dtype = 'float')        # rotation matrix
        C = np.array([[0,0,self.M[2,2]*v],[0,0,self.M[1,1]*u],[self.M[2,2]*v,-self.M[1,1]*u,0]],dtype = 'float')
        F = np.dot(np.linalg.inv(self.M),(-np.dot(self.D,nu)-np.dot(C,nu)+tw))
        fu = F[0]
        fr = F[2]
        
        nu_dot = F+np.dot(np.linalg.inv(self.M),np.dot(self.B,Tf))
        eta_dot = np.dot(R,nu)+vc
        
        nu = ts*nu_dot+nu
        eta = ts*eta_dot+eta
        
        return nu,eta,nu_dot,eta_dot,Tf,fu,fr
    
    def state_update2(self,nu,eta,Tf,Tc,vc,tw,ts):
        '''
        

        Parameters
        ----------
        nu : TYPE
            velocity state of USV
        eta : TYPE
            position state of USV
        Tf : TYPE
            thrust filtered value
        Tc : TYPE
            thrust command input
        vc : TYPE
            currents
        tw : TYPE
            unknown disturbances
        ts : TYPE
            simulation step

        Returns
        -------
        next_nu
        next_eta
        nu_dot
        eta_dot
        next_Tf
        fu
        fr

        '''
        Tf_dot = -(Tf-Tc)/self.Tf_time
        
        # limit the change rate of thrust command
        for i in range(2):
            if Tf_dot[i] >= self.Tdot_max:
                Tf_dot[i] = self.Tdot_max
            if Tf_dot[i] <= self.Tdot_min:
                Tf_dot[i] = self.Tdot_min
                
        # update the thrust
        Tf = ts*Tf_dot+Tf
        
        # limit the amplitude of thrust command
        for i in range(2):
            if Tf[i] >= self.T_max:
                Tf[i] = self.T_max
            if Tf[i] <= self.T_min:
                Tf[i] = self.T_min
                
        # USV's state update
        psi = eta[2]
        u = nu[0]
        v = nu[1]
        R = np.array([[np.cos(psi),-np.sin(psi),0],[np.sin(psi),np.cos(psi),0],[0,0,1]],dtype = 'float')        # rotation matrix
        C = np.array([[0,0,self.M[2,2]*v],[0,0,self.M[1,1]*u],[self.M[2,2]*v,-self.M[1,1]*u,0]],dtype = 'float')
        F = np.dot(np.linalg.inv(self.M),(-np.dot(self.D,nu)-np.dot(C,nu)+tw))
        fu = F[0]
        fv = F[1]
        fr = F[2]
        tau = np.dot(self.B,Tf)
        U = np.sqrt(u**2+v**2)
        if U==0:
            U=0.1
        fU = (u*fu+v*fv+(u-U)*tau[0]/self.M[0,0])/U
        
        nu_dot = F+np.dot(np.linalg.inv(self.M),np.dot(self.B,Tf))
        eta_dot = np.dot(R,nu)+vc
        
        nu = ts*nu_dot+nu
        eta = ts*eta_dot+eta
        
        return nu,eta,nu_dot,eta_dot,Tf,fU,fr

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    