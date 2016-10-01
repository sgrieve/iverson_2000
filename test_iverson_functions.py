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
    
    t_stars = np.linspace(0.1, 10000, 10)
    
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
    
    R2 = IvF.R_fn(0.0001)
    print "R2 is"
    print R2    
  
    # These are some characteristics of the slope
    alpha = math.radians(15.)
    
    # Here are the rainfall intensities
    Iz_over_Kz_steady = 0.1     # the intensity of the "steady state" pressure profile
    Iz_over_Kz = 1              # intensity of storm event
    
    # set up the spatial coordinates
    Z = np.linspace(0.0, 6., 13)
    beta = IvF.Beta_fn(alpha, Iz_over_Kz_steady)   
    
    #==========================================================
    # test the part of the function that calculates pressure
    print "beta is: "
    print beta
    
    t_star = 0.0001
    T_star = 10
    
    # This is the depth of the steady state water table
    d = 2
    
    Z = 2.3    
    
    # Test the pressure function. This solves equations 27a and b
    this_psi = IvF.psi(Z,beta,d,Iz_over_Kz,t_star,T_star)
    print "This Z is: "
    print Z
    print "This psi is: "
    print this_psi
    #===========================================================
 
    # time and peak time of rainfall duration are in weeks so we need to 
    # convert them to seconds 
    t= 2
    T = 10
    t_sec = IvF.weeks_to_secs(t)
    T_sec = IvF.weeks_to_secs(T)
    Do = 0.000001
    D_hat = IvF.D_hat_fn(Do, alpha)   
    
    Zs = np.linspace(0.01, 6., 10)
    print "Z is: "
    print Zs
    print "Now hold on a sec whi I calculate pressure"    
    this_psi = IvF.psi_dimensional_t(Zs,beta,d,Iz_over_Kz,D_hat,t,T)

    
    print "And now for the dimensional psi: "
    print this_psi
    
    
    # reset the t variable for the figure
    t = [0,4,8,12,24]
 
    # calculate figure 7
    print "I am making figure 7b"
    IvF.Iverson_Fig_7(t, T, d, Do, alpha, Iz_over_Kz, Iz_over_Kz_steady)
    
    print "Okay now I am going to do figure 7a"
    # now do the other one
    fig7a_ts = [0,2,6,10,20]
    t = IvF.days_to_weeks(fig7a_ts)
    fig_7a_T = 10
    T = IvF.days_to_weeks(fig_7a_T)
    IvF.Iverson_Fig_7(t, T, d, Do, alpha, Iz_over_Kz, Iz_over_Kz_steady)
    
if __name__ == "__main__":
    compare_linear_to_loop()  
    