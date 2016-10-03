import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 18
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'


def days_to_secs(days):
    '''
    Convert a number of days to a number of seconds.
    '''
    return days * 24. * 60. * 60.


def weeks_to_secs(weeks):
    '''
    Convert a number of weeks to a number of seconds.
    '''
    return days_to_secs(weeks * 7.)


def erfc(XX):
    '''
    Apply the math.erfc() function to an array or list of values, returning a
    numpy array.

    Not currently used, but may be useful.
    '''
    X = np.asarray(XX)
    data = []
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
    return (x * np.sin(alpha)) + (z * np.cos(alpha))


def D_hat_fn(Do, alpha):
    '''
    Equation 26c - compute an effective hydraulic diffusivity.

    Do is a reference diffusivity and alpha is the slope gradient.
    '''
    return 4. * Do * (np.cos(alpha) ** 2.)


def Iz_over_Kz_steady(alpha, Iz_over_Kz):
    '''
    Calculate the Iz_over_Kz_steady constant used in beta.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8, but he also defines
    it as a function of Iz_over_Kz near the top of the first column on page 1902
    '''
    return np.cos(alpha) * Iz_over_Kz


def Beta_fn(alpha, Iz_over_Kz_steady):
    '''
    Calculate the beta constant used in equation 27a and b.

    Alpha is the gradient and Iz_over_Kz_steady is Steady state vertical water
    influx rate, Iverson sets this to 0.1 in Figures 7 and 8.
    '''
    return (np.cos(alpha) ** 2.) - Iz_over_Kz_steady


def t_T_star_fn(t, D_hat, Z):
    '''
    Equation 27c and 27d to nondimensionalise the time parameters.

    D_hat is an effective hydraulic diffusivity and Z is the elevation head.

    Both have the same form, but 27c takes t (time) whereas d takes T (rainfall
    duration).
    '''
    return t / ((Z * Z) * D_hat)


def R_fn(t_star):
    '''
    Compute the response function for a t* value.

    Equation 27e
    '''

    return (np.sqrt(t_star / np.pi) * np.exp(-1. / t_star) -
            math.erfc(1. / np.sqrt(t_star)))


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


def Iverson_Eq_27ab(t, T, Do, alpha, Z, Iz_over_Kz=1., Iz_over_Kz_steady=0.1):
    '''
    Equations 27 a and b

    Does not work.
    '''

    d = 2.

    D_hat = D_hat_fn(Do, alpha)
    t_star = t_T_star_fn(t, D_hat, Z)

    T_star = t_T_star_fn(T, D_hat, Z)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)

    if t <= T:
        return ((beta * (1. - (d / Z))) + (Iz_over_Kz * R_fn(t_star)))
    else:
        return ((beta * (1. - (d / Z))) + (Iz_over_Kz * (R_fn(t_star) -
                R_fn(t_star - T_star))))


def Iverson_Fig_7(t, T, Do, alpha, Iz_over_Kz, Iz_over_Kz_steady):
    '''
    Reproduces Figure 7. Currently does not work.
    '''

    Zs = np.linspace(0.01, 6., 10)
    beta = Beta_fn(alpha, Iz_over_Kz_steady)

    psi = []
    beta_line = []

    for Z in Zs:
        psi_ = Iverson_Eq_27ab(t, T, Do, alpha, Z)
        print Z, psi_ * Z
        psi.append(psi_ * Z)
        beta_line.append(beta * Z)

    print "beta_line is"
    print beta_line

    plt.gca().invert_yaxis()
    plt.plot(psi, Zs)
    plt.plot(beta_line, Zs, 'k--', label='$\\beta$ Line')
    # plt.xlim(-2, 5)
    legend = plt.legend()
    legend.get_frame().set_linewidth(0.)
    plt.xlabel('Pressure head (m)')
    plt.ylabel('Depth (m)')
    plt.tight_layout()
    plt.show()


#Iverson_Fig_5(0.1)
#Iverson_Fig_7(6., 10., 0.000001, math.radians(15.), 1., 0.1)
# Iverson_Fig_6()

#a = [1,2,3,4,5]
#t = np.asarray(a)
#R = R_fn(t)
#print a
