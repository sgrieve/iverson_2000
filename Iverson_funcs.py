import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams


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
rcParams['legend.labelspacing'] =0.1
rcParams['legend.columnspacing'] =0.1 


def days_to_secs(days):
    '''
    Convert a number of days to a number of seconds.
    '''
    days = np.asarray(days)
    if days.size == 1:
        return days * 24. * 60. * 60.
    else:
        return np.multiply(days,86400.)


def weeks_to_secs(weeks):
    '''
    Convert a number of weeks to a number of seconds.
    '''
    weeks = np.asarray(weeks)
    if weeks.size == 1:
        return days_to_secs(weeks * 7.)
    else:
        return np.multiply(days_to_secs(weeks),7.)


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


def Z_fn(x, z, alpha):
    '''
    Calculate Z, the elevation head, from the equations given in Figure 3's
    caption.

    Alpha is the gradient, x and z are coordinates.

    If Z and z share an orign, it simplifies to:
    return z / np.cos(alpha)
    '''
    term1 = np.multiply(x,np.sin(alpha))
    term2 = np.multiply(z,np.cos(alpha))
    return ( np.add(term1,term2))


def D_hat_fn(Do, alpha):
    '''
    Equation 26c - compute an effective hydraulic diffusivity.

    Do is a reference diffusivity and alpha is the slope gradient.
    '''
    return 4. * Do * (np.cos(alpha) ** 2.)

def Iz_over_Kz_steady(alpha,Iz_over_Kz):
    '''
    Calculate the Iz_over_Kz_steady constant used in beta.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8, but he also defines
    it as a function of Iz_over_Kz near the top of the first column on page 1902
    ''' 
    return np.cos(alpha)*Iz_over_Kz

def Beta_fn(alpha, Iz_over_Kz_steady):
    '''
    Calculate the beta constant used in equation 27a and b.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8.
    '''
    return (np.cos(alpha) ** 2.) - Iz_over_Kz_steady

def Beta_line(Z,beta):
    '''
    Returns the "beta line", which is the maximum pore pressure
    '''
    beta_line = beta*Z
    return beta_line
        


def t_T_star_fn(t, D_hat, Z):
    '''
    Equation 27c and 27d to nondimensionalise the time parameters.

    D_hat is an effective hydraulic diffusivity and Z is the elevation head.

    Both have the same form, but 27c takes t (time) whereas d takes T (rainfall
    duration).
    '''
    Z2 = np.multiply(Z,Z)
    denominator = np.multiply(Z2,D_hat)
          
    return np.divide(t,denominator)


def R_fn(t_star):
    '''
    Compute the response function for a t* value.

    Equation 27e
    '''
    
    multiple_bit = np.multiply(np.sqrt(t_star / np.pi),np.exp(-1. / t_star)) 
    one_ov_sqrt_tstar = 1. / np.sqrt(t_star)
    
    #print "Numbers are (mulitply, one ov tstar): "
    #print multiple_bit
    #print one_ov_sqrt_tstar
    
    #print "is this the problem??"
    #print array_erfc(one_ov_sqrt_tstar)
    
    #print "R_fn is"    
    #print (np.subtract(multiple_bit,array_erfc(one_ov_sqrt_tstar)))
    return np.subtract(multiple_bit,array_erfc(one_ov_sqrt_tstar))

def psi(Z,beta,d,Iz_over_Kz,t_star,T_star):
    '''
    Compute psi from equation 27a and b
    This is slightly confusing because the dimensional time is nondimensionalised
    By depth (see equations 27c,d), so that for a single t_star or T_star value,
    the dimensional time varies with depth
    '''
    
    # This is effectively the steady state water table (or initial condition)
    # Iverson seems to just pull this out of the air. Perhaps it is based on
    # measurements?
    first_term = beta*(Z-d)
    #print "first_term is: "
    #print first_term

    #print "For a t_star of: "+str(t_star)+" R_fn is"
    #print R_fn(t_star)    
    
    # This solves the equation, based on the response function (R_fn), 
    # which is equation 27e
    if t_star < T_star:
        second_term = Z*Iz_over_Kz*R_fn(t_star)
        #print "Second term is: " + str(second_term)
    else:
        second_term = Z*Iz_over_Kz*(R_fn(t_star)-R_fn(t_star-T_star))
        #print "Second term is (t_star > T_star): " + str(second_term)
        
    psi = first_term+second_term
    
    return psi
            

    
def psi_dimensional_t(Z,beta,d,Iz_over_Kz,D_hat,t,T):
    '''
    Compute psi from equation 27a and b, but using dimensional time
    A bit slow since I haven't vectorised the calculations. 
    '''

    psi_dimensional = []
    for z in Z:
        # first get the nondimensional time. Note that according to 
        # equations 27c,d the dimensionless time is a function of depth, 
        # so each point below the surface has a different t_star and T_star
        t_star = t*D_hat/(z*z)
        T_star = T*D_hat/(z*z)
        
        # Now calculate psi on the basis of the dimensional psi
        this_psi = psi(z,beta,d,Iz_over_Kz,t_star,T_star)
        psi_dimensional.append(this_psi)

                           
    return psi_dimensional
    


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

def Iverson_Fig_7(t, T, d, Do, alpha, Iz_over_Kz, Iz_over_Kz_steady):
    '''
    Reproduces Figure 7. Currently does not work.
    Supply t as a vector in weeks and T in weeks
    '''

    # Get parameters for psi curves
    D_hat = D_hat_fn(Do, alpha)   
    Zs = np.linspace(0.01, 6., 10)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)
    beta_line = Beta_line(Zs,beta)
    
    initial_line = psi_dimensional_t(Zs,beta,d,Iz_over_Kz,D_hat,0,T)    

    # times in weeks
    ts_in_weeks = [0,4,8,12,24]
    ts = weeks_to_secs(ts_in_weeks)

    Fig1 = plt.figure(1, facecolor='white',figsize=(3.26,3.26))  

    Fig1.gca().invert_yaxis()   
    plt.plot(beta_line, Zs, 'k--', label='$\\beta$ Line')

    Ts = weeks_to_secs(T)
    
    for t_week in t:
        ts = weeks_to_secs(t_week)
        this_label = 't = '+str(t_week)+' weeks'
        psi = psi_dimensional_t(Zs,beta,d,Iz_over_Kz,D_hat,ts,Ts)    
        plt.plot(psi, Zs, label=this_label)

    # plt.xlim(-2, 5)
    legend = plt.legend()
    legend.get_frame().set_linewidth(0.)
    plt.xlabel('Pressure head (m)')
    plt.ylabel('Depth (m)')
    plt.title('T = '+str(T)+' weeks')
    plt.tight_layout()
    plt.show()



