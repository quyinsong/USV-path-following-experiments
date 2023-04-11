# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 15:13:07 2022

@author: qys
"""

def euler2(xdot,x,ts):
    x = x+xdot*ts
    return x