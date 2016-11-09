from __future__ import print_function
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patches as patches


label_size = 8
axis_size = 12

# Set up fonts for plots
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size
rcParams['xtick.major.size'] = 4
rcParams['ytick.major.size'] = 4
rcParams['legend.fontsize'] = label_size
rcParams['legend.handletextpad'] = 0.05
rcParams['legend.labelspacing'] = 0.1
rcParams['legend.columnspacing'] = 0.1


#==============================================================================
def minutes_to_secs(minutes):
    '''
    Convert a number of days to a number of seconds.
    '''
    minutes = np.asarray(minutes)
    if minutes.size == 1:
        return minutes * 60.
    else:
        return np.multiply(minutes, 60.)
#==============================================================================



#==============================================================================
def days_to_secs(days):
    '''
    Convert a number of days to a number of seconds.
    '''
    days = np.asarray(days)
    if days.size == 1:
        return days * 24. * 60. * 60.
    else:
        return np.multiply(days, 86400.)
#==============================================================================
        
        
#==============================================================================
def weeks_to_secs(weeks):
    '''
    Convert a number of weeks to a number of seconds.
    '''
    weeks = np.asarray(weeks)
    if weeks.size == 1:
        return days_to_secs(weeks * 7.)
    else:
        return np.multiply(days_to_secs(weeks), 7.)
#==============================================================================
        
        
#==============================================================================        
def days_to_weeks(days):
    '''
    Convert a number of days to a number of weeks.
    '''
    days = np.asarray(days)
    if days.size == 1:
        return days / 7.
    else:
        return np.divide(days, 7.)
#==============================================================================
        
        
        
#==============================================================================
def array_erfc(XX):
    '''
    Apply the math.erfc() function to an array or list of values, returning a
    numpy array.

    Not currently used, but may be useful.
    '''

    X = np.asarray(XX)
    #print X.size
    #print "The array is:"
    #print X

    data = []
    if X.size == 1:
        return math.erfc(X)
    else:
        for x in X:
            data.append(math.erfc(x))
            return np.array(data)
#==============================================================================
            
            
            
#==============================================================================
def Z_fn(x, z, alpha):
    '''
    Calculate Z, the elevation head, from the equations given in Figure 3's
    caption.

    Alpha is the gradient, x and z are coordinates.

    If Z and z share an orign, it simplifies to:
    return z / np.cos(alpha)
    '''
    term1 = np.multiply(x, np.sin(alpha))
    term2 = np.multiply(z, np.cos(alpha))
    return (np.add(term1, term2))
#==============================================================================
    
    
    
#==============================================================================
def D_hat_fn(Do, alpha):
    '''
    Equation 26c - compute an effective hydraulic diffusivity.

    Do is a reference diffusivity and alpha is the slope gradient.
    '''
    return 4. * Do * (np.cos(alpha) ** 2.)
#==============================================================================
    
    
    
#==============================================================================
def Iz_over_Kz_steady(alpha, Iz_over_Kz):
    '''
    Calculate the Iz_over_Kz_steady constant used in beta.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8, but he also defines
    it as a function of Iz_over_Kz near the top of the first column on page 1902
    '''
    return np.cos(alpha) * Iz_over_Kz
#==============================================================================
    
    
    
#==============================================================================
def Beta_fn(alpha, Iz_over_Kz_steady):
    '''
    Calculate the beta constant used in equation 27a and b.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8.
    '''
    return (np.cos(alpha) ** 2.) - Iz_over_Kz_steady
#==============================================================================
    
    
    
#==============================================================================
def Beta_line(Z, beta):
    '''
    Returns the "beta line", which is the maximum pore pressure
    '''
    beta_line = beta * Z
    return beta_line
#==============================================================================
    
def Correct_psi(Z,Psi,beta):
    '''
    Reduces Psi to the beta line, which is the theoretical maximum pore pressure
    '''   

    #print("Psi values are: ")
    #print(Psi)    
    
    for i,Psi_val in enumerate(Psi):
        this_beta = Z[i]*beta
        if Psi_val > this_beta:
            Psi[i] = this_beta
            
    return Psi
    
    
    
#==============================================================================
def t_T_star_fn(t, D_hat, Z):
    '''
    Equation 27c and 27d to nondimensionalise the time parameters.

    D_hat is an effective hydraulic diffusivity and Z is the elevation head.

    Both have the same form, but 27c takes t (time) whereas d takes T (rainfall
    duration).
    '''
    Z2 = np.multiply(Z, Z)
    denominator = np.multiply(Z2, D_hat)

    return np.divide(t, denominator)
#==============================================================================
    
    
    
#==============================================================================
def R_fn(t_star):
    '''
    Compute the response function for a t* value.

    Equation 27e
    '''

    multiple_bit = np.multiply(np.sqrt(t_star / np.pi), np.exp(-1. / t_star))
    one_ov_sqrt_tstar = 1. / np.sqrt(t_star)

    #print "Numbers are (mulitply, one ov tstar): "
    #print multiple_bit
    #print one_ov_sqrt_tstar

    #print "is this the problem??"
    #print array_erfc(one_ov_sqrt_tstar)

    #print "R_fn is"
    #print (np.subtract(multiple_bit,array_erfc(one_ov_sqrt_tstar)))
    return np.subtract(multiple_bit, array_erfc(one_ov_sqrt_tstar))
#==============================================================================
    
    
    
#==============================================================================
def psi(Z, beta, d, Iz_over_Kz, t_star, T_star):
    '''
    Compute psi from equation 27a and b
    This is slightly confusing because the dimensional time is nondimensionalised
    By depth (see equations 27c,d), so that for a single t_star or T_star value,
    the dimensional time varies with depth
    '''

    # This is effectively the steady state water table (or initial condition)
    # Iverson seems to just pull this out of the air. Perhaps it is based on
    # measurements?
    first_term = beta * (Z - d)
    #print "first_term is: "
    #print first_term

    #print "For a t_star of: "+str(t_star)+" R_fn is"
    #print R_fn(t_star)

    # This solves the equation, based on the response function (R_fn),
    # which is equation 27e
    if t_star < T_star:
        second_term = Z * Iz_over_Kz * R_fn(t_star)
        #print "Second term is: " + str(second_term)
    else:
        second_term = Z * Iz_over_Kz * (R_fn(t_star) - R_fn(t_star - T_star))
        #print "Second term is (t_star > T_star): " + str(second_term)

    psi = first_term + second_term

    return psi
#==============================================================================
    
#==============================================================================
def psi_transient(Z, beta, d, Iz_over_Kz, t_star, T_star):
    '''
    Compute psi from equation 27a and b
    This one only has the transient component
    '''

    # This solves the equation, based on the response function (R_fn),
    # which is equation 27e
    if t_star < T_star:
        second_term = Z * Iz_over_Kz * R_fn(t_star)
        #print "Second term is: " + str(second_term)
    else:
        second_term = Z * Iz_over_Kz * (R_fn(t_star) - R_fn(t_star - T_star))
        #print "Second term is (t_star > T_star): " + str(second_term)

    psi = second_term

    return psi
#============================================================================== 

       
#==============================================================================
def psi_dimensional_t(Z, beta, d, Iz_over_Kz, D_hat, t, T):
    '''
    Compute psi from equation 27a and b, but using dimensional time
    A bit slow since I haven't vectorised the calculations.
    Times need to be in seconds
    '''

    psi_dimensional = []
    for z in Z:
        # first get the nondimensional time. Note that according to
        # equations 27c,d the dimensionless time is a function of depth,
        # so each point below the surface has a different t_star and T_star
        t_star = t * D_hat / (z * z)
        T_star = T * D_hat / (z * z)

        # Now calculate psi on the basis of the dimensional psi
        this_psi = psi(z, beta, d, Iz_over_Kz, t_star, T_star)
        psi_dimensional.append(this_psi)

    return psi_dimensional
#==============================================================================
    

#==============================================================================
def psi_dimensional_t_transient(Z, beta, d, Iz_over_Kz, D_hat, t, T):
    '''
    Compute psi from equation 27a and b, but using dimensional time
    A bit slow since I haven't vectorised the calculations.
    Only calculates the transient component of psi for use with 
    time series of rainfall
    times need to be in seconds
    '''

    psi_dimensional = []
    for z in Z:
        # first get the nondimensional time. Note that according to
        # equations 27c,d the dimensionless time is a function of depth,
        # so each point below the surface has a different t_star and T_star
        t_star = t * D_hat / (z * z)
        T_star = T * D_hat / (z * z)

        # Now calculate psi on the basis of the dimensional psi
        this_psi = psi_transient(z, beta, d, Iz_over_Kz, t_star, T_star)
        psi_dimensional.append(this_psi)

    return psi_dimensional
#==============================================================================

#==============================================================================
# This loads a time series that has rain two columns, a rainfall duration
# and a rainfall rate. 
# It then superposes these to get the Psi value at a given time
#==============================================================================
def psi_dimensional_from_time_series(durations,intensities,Z, beta, d, D_hat, t):
    
    # this is the steady component of psi. Basically transient pressure builds 
    # on top of this. Not clear where Iverson comes up with these numbers. 
    steady_psi = beta * (Z - d)
    cumulative_psi = np.asarray(steady_psi)
        
    
    # Now we try to construct the transient pressure. 
    # loop through the record getting cumulative times
    starting_times = []
    starting_times.append(0)
    cumulative_time = 0
    count = 0
    end_count_found = False
    end_count = 0
    for this_duration in durations:
        cumulative_time = cumulative_time+this_duration
        
        # the cumulative time is the time at the end of this timestep. 
        # if the cumulative time  is less than the time of simulation, 
        # then we need to acount for this pulse of rainfall        
        if t < cumulative_time:
            if end_count_found == False:
                end_count_found = True
                end_count = count
        
        count = count+1
        starting_times.append(cumulative_time)
    
    # we don't need the last element
    del starting_times[-1]
    
    # If we didn't find the end count it means the rainfall records have ended and we need
    # all of the data        
    if end_count_found == False:
        end_count = count-1     # The minus one is needed since we have counted past the end of the index
        

    # okay, now get the transients from superposition 
    # First we need to figure out how many of these we will need
    print("end count is: " + str(end_count))
    
    for i,time in enumerate(starting_times):
        print("time is: "+str(time))
        
        eff_t = t-time
        this_intensity = intensities[i]
        this_duration = durations[i]
        
        print("Eff T: "+str(eff_t)+" and intensity is: "+str(this_intensity)+" and duration is: " +str(this_duration))
        
        # now get the psi values.
        this_psi = psi_dimensional_t_transient(Z, beta, d, this_intensity, D_hat, eff_t, this_duration)
        
        cumulative_psi = np.add(cumulative_psi,this_psi)
        
        print("this psi is:")
        print(this_psi)
        
        print ("cumulative psi is: ")
        print(cumulative_psi)



    
    
#==============================================================================
def F_f(alpha, friction_angle):
    '''
    Compute F_f from equation 28
    '''

    tan_alpha = np.tan(alpha)
    tan_friction_angle = np.tan(friction_angle)

    return np.divide(tan_friction_angle, tan_alpha)
#==============================================================================
    
    
#==============================================================================    
def F_c(cohesion, weight_of_soil, Z, alpha):
    '''
    Compute F_c from equation 28
    '''   
    denominator_1 = np.multiply(weight_of_soil, Z)
    denominator_2 = np.multiply(np.sin(alpha), np.cos(alpha))
    denominator = np.multiply(denominator_1, denominator_2)
    
    #print("Denominator 1 is:")
    #print(denominator_1)
    #print("Denominator 2 is:")
    #print(denominator_2)
    #print("Entire denominator: ")
    #print(denominator)

    return np.divide(cohesion, denominator)
#==============================================================================
    
    
#==============================================================================
def F_w(psi_val, weight_of_water, weight_of_soil, alpha, friction_angle, Z):
    '''
    Compute F_w from equation 28
    The Psi value is the dimensional Psi at a dimensional time. 
    '''  
    numerator_1 = np.multiply(psi_val, weight_of_water)
    numerator_2 = np.multiply(numerator_1, np.tan(friction_angle))
    numerator = np.multiply(-1., numerator_2)

    denominator_1 = np.multiply(weight_of_soil, Z)
    denominator_2 = np.multiply(np.sin(alpha), np.cos(alpha))
    denominator = np.multiply(denominator_1, denominator_2)

    return np.divide(numerator, denominator)
#==============================================================================


#==============================================================================   
def FS(psi_val, weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, Z):
    '''
    Compute FS from equation 28a (same as equation 29, but broken into different components)
    The Psi value is the dimensional Psi at a dimensional time. 
    '''
    this_F_f = F_f(alpha, friction_angle)
    this_F_c = F_c(cohesion, weight_of_soil, Z, alpha)
    FS0 = np.add(this_F_f,this_F_c)
    FSprime = F_w(psi_val, weight_of_water, weight_of_soil, alpha, friction_angle, Z)
    
    FS = np.add(FS0,FSprime)
    return FS
#============================================================================== 
    
    
#============================================================================== 
def FS_fxn_t_T_Z(Zs, t_sec, T_sec, weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, beta, d, Iz_over_Kz, D_0):
    '''
    Compute FS from equation 28a (same as equation 29, but broken into different components)
    The Psi value is the dimensional Psi at a dimensional time. 
    '''   
    # Get the normalised diffusivity
    D_hat = D_hat_fn(D_0, alpha)    
    
    # get pressure
    this_psi = psi_dimensional_t(Zs, beta, d, Iz_over_Kz, D_hat, t_sec, T_sec)

    # Correct Psi: it is limited by the beta curve (which is just the saturated pore pressure)    
    # NOTE IVERSON DOESN'T USE THIS IN HIS FIGURES EVEN THOUGH HE SAYS YOU SHOULD    
    corr_psi = Correct_psi(Zs,this_psi,beta)
    
    # This is what iverson does to generate figres 10 and 11
    #corr_psi = this_psi
    
    # Now get the Factor of safety
    this_FS = FS(corr_psi, weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, Zs)
    
    return this_FS
#============================================================================== 





   

#==============================================================================
def Iverson_Fig_5(T_star):
    '''
    Reproduces Figure 5. Pass in T* values of 0.1, 1.0 and 10 to generate
    subplots A, B and C.
    '''

    t_stars = np.linspace(0.1, 1000, 10000)

    vals = []

    for t in t_stars:
        if t <= T_star:
            vals.append(R_fn(t))
        else:
            vals.append(R_fn(t) - R_fn(t - T_star))

    plt.plot(t_stars, vals)

    plt.vlines(T_star, 0., max(vals), linestyle='--',
               label='$T^*$ = {0}'.format(T_star))

    ax = plt.gca()
    ax.set_xscale("log", nonposx='clip')
    plt.xlabel('Normalised Time, $t^*$')
    plt.ylabel('R')
    legend = plt.legend()
    legend.get_frame().set_linewidth(0.)
    plt.tight_layout()
    plt.savefig('Fig_5.png')
#==============================================================================
    



    
    
#==============================================================================
def Iverson_Fig_6():
    '''
    Reproduces Figure 6.
    Calls the response function for a range of T* values and identifies the
    peak R value and plots the peak t* and peak R values.
    '''

    T_Stars = np.linspace(0.001, 1000, 1000)
    t_stars = np.linspace(0.01, 1000, 1000)

    t_peak = []
    R_peak = []

    for T in T_Stars:
        vals = []

        for t in t_stars:
            if t <= T:
                vals.append(R_fn(t))
            else:
                vals.append(R_fn(t) - R_fn(t - T))

        index = np.argmax(vals)

        t_peak.append(t_stars[index])
        R_peak.append(vals[index])

    plt.plot(T_Stars, t_peak, 'k--', label='$t^*$ peak')
    plt.plot(T_Stars, R_peak, 'r-', label='R peak')

    ax = plt.gca()
    ax.set_xscale("log", nonposx='clip')
    ax.set_yscale("log", nonposy='clip')
    plt.xlabel('Normalised Rainfall Duration, $T^*$')
    plt.ylabel('Normalised Peak Time and Peak Response')
    plt.ylim(0.001, 1000.)
    legend = plt.legend(loc=2)
    legend.get_frame().set_linewidth(0.)
    plt.tight_layout()
    plt.savefig('Fig_6.png')
#==============================================================================
    
    
    
#==============================================================================
def Iverson_Fig_7(t, T, d, Do, alpha, Iz_over_Kz, Iz_over_Kz_steady):
    '''
    Reproduces Figure 7. Currently does not work.
    Supply t as a vector in weeks and T in weeks
    '''

    # Get parameters for psi curves
    D_hat = D_hat_fn(Do, alpha)
    Zs = np.linspace(0.01, 6., 10)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)
    beta_line = Beta_line(Zs, beta)

    initial_line = psi_dimensional_t(Zs, beta, d, Iz_over_Kz, D_hat, 0, T)

    # times in weeks
    ts_in_weeks = [0, 4, 8, 12, 24]
    ts = weeks_to_secs(ts_in_weeks)

    Fig1 = plt.figure(1, facecolor='white', figsize=(10, 8))

    Fig1.gca().invert_yaxis()
    plt.plot(beta_line, Zs, 'k--', label='$\\beta$ Line')

    Ts = weeks_to_secs(T)

    for t_week in t:
        ts = weeks_to_secs(t_week)
        this_label = 't = ' + str(t_week) + ' weeks'
        psi = psi_dimensional_t(Zs, beta, d, Iz_over_Kz, D_hat, ts, Ts)
        plt.plot(psi, Zs, label=this_label)

    # plt.xlim(-2, 5)
    legend = plt.legend()
    legend.get_frame().set_linewidth(0.)
    plt.xlabel('Pressure head (m)')
    plt.ylabel('Depth (m)')
    plt.title('T = ' + str(T) + ' weeks')
    plt.tight_layout()
    plt.savefig("IversonFig7_" + str(int(T)) + ".png", format="png")
    plt.show()
#==============================================================================
    
    
#==============================================================================
def Iverson_FoS_Fig10(weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, d, Iz_over_Kz, Iz_over_Kz_steady, D_0, t_sec, T_sec,name_string):
    
    # Get parameters for psi curves
    Zs = np.linspace(0.001, 6., 200)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)

    Fig1 = plt.figure(1, facecolor='white', figsize=(10, 8))

    Fig1.gca().invert_yaxis()
   

    for this_time in t_sec:
        
        this_FS = FS_fxn_t_T_Z(Zs, this_time, T_sec, weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, beta, d, Iz_over_Kz, D_0)
        
        this_label = 't = ' + str(this_time) + ' seconds'        
        plt.plot(this_FS, Zs, label=this_label,linewidth=3)
    
    ax1 = plt.gca()    
    ax1.add_patch(patches.Rectangle(
        (0.9, 0),   # (x,y)
        0.1,          # width
        6,          # height
        alpha=0.2,
        facecolor="red",
        linewidth=2,
    ))
    plt.text(0.92, 2, 'Unstable', rotation=90, fontsize=24)
    
    #ax = plt.gca()
    plt.xlim(0.9, 1.4)
    
    legend = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=2.,fontsize = 20)        
    legend.get_frame().set_linewidth(1.)
    plt.xlabel('Factor of Safety (dimensionless)')
    plt.ylabel('Depth (m)')
    #plt.title('T = ' + str(T) + ' weeks')
    plt.tight_layout()    
    plt.savefig(name_string, format="png")
    plt.show()       
#==============================================================================    
 
#==============================================================================
def Iverson_FoS_Fig11(weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, d, Iz_over_Kz, Iz_over_Kz_steady, D_0, t_sec, T_sec,name_string):
    
    # Get parameters for psi curves
    Zs = np.linspace(0.001, 0.7, 200)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)

    Fig1 = plt.figure(1, facecolor='white', figsize=(10, 8))

    Fig1.gca().invert_yaxis()
   

    for this_time in t_sec:
        
        this_FS = FS_fxn_t_T_Z(Zs, this_time, T_sec, weight_of_water, weight_of_soil, alpha, cohesion, friction_angle, beta, d, Iz_over_Kz, D_0)
        
        this_label = 't = ' + str(this_time) + ' seconds'        
        plt.plot(this_FS, Zs, label=this_label,linewidth=3)
    
    ax1 = plt.gca()    
    ax1.add_patch(patches.Rectangle(
        (0.5, 0),   # (x,y)
        0.5,          # width
        0.7,          # height
        alpha=0.2,
        facecolor="red",
        linewidth=2,
    ))
    plt.text(0.92, 2, 'Unstable', rotation=90, fontsize=24)
    
    #ax = plt.gca()
    plt.xlim(0.5, 2.0)
    
    legend = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=2.,fontsize = 20)        
    legend.get_frame().set_linewidth(1.)
    plt.xlabel('Factor of Safety (dimensionless)')
    plt.ylabel('Depth (m)')
    #plt.title('T = ' + str(T) + ' weeks')
    plt.tight_layout()    
    plt.savefig(name_string, format="png")
    plt.show()       
#==============================================================================    
   