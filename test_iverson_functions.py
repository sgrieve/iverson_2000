# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 09:22:28 2016

@author: smudd
"""

import Iverson_funcs as IvF
import iverson_2000 as I2000
import numpy as np
import math

def compare_linear_to_loop():
    
    t_stars = np.linspace(0.1, 1000, 10)
    
    print "t_stars are:"
    print t_stars

    vals = []
    vals2 = []

    for t in t_stars:
        vals.append(I2000.R_fn(t))
        vals2.append(IvF.R_fn(t))
    

    R = IvF.R_fn(t_stars)
    
    print "Version 1: "
    print vals
    
    print "Version 2: "
    print R
    
    print "Version 3: "
    print vals2
    
    R2 = IvF.R_fn(0.1)
    print "R2 is"
    print R2
    
    
  
    alpha = math.radians(15.)
    Iz_over_Kz_steady = 0.1
    
    Z = np.linspace(0.5, 6., 12)
    beta = IvF.Beta_fn(alpha, Iz_over_Kz_steady)   
    
    print "beta is: "
    print beta
    
    t_star = 0.00000000001
    T_star = 10
    Iz_over_Kz = 1
    d = 2
    
    this_psi = IvF.psi(Z,beta,d,Iz_over_Kz,t_star,T_star)
    print "This Z is: "
    print Z
    print "This psi is: "
    print this_psi
  
    
if __name__ == "__main__":
    compare_linear_to_loop()  
    