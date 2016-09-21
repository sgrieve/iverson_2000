import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rcParams


rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = 18
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'


def erfc(x):
    data = []
    for x in X:
        data.append(math.erfc(x))

    return np.array(data)


def R_fn(t_star):
    return (np.sqrt(t_star / np.pi) * np.exp(-1. / t_star) -
            math.erfc(1. / np.sqrt(t_star)))


def Beta_fn(alpha, Iz_over_Kz):
    return (np.cos(alpha) ** 2.) * Iz_over_Kz


def Z_fn(x, z, alpha):
    return (x * np.sin(alpha)) + (z * np.cos(alpha))

    # if Z and z share an orign, the above simplifies to:
    # return z / np.cos(alpha)


def D_hat_fn(Do, alpha):
    return 4. * Do * (np.cos(alpha) ** 2.)


def t_T_star_fn(t, D_hat, Z):
    return t / ((Z * Z) * D_hat)


def days_to_secs(days):
    return days * 24. * 60. * 60.


def weeks_to_secs(weeks):
    return days_to_secs(weeks * 7.)


def Iverson_Eq_27ab(t, T, Do, alpha, d, Iz_over_Kz=1.):
    Z = 6.
    x = 10.
    #z = d
    #Z = d
    #d = 2.
    #Z = Z_fn(x, z, alpha)

    D_hat = D_hat_fn(Do, alpha)
    t_star = t_T_star_fn(t, D_hat, Z)
    T_star = t_T_star_fn(T, D_hat, Z)
    beta = Beta_fn(alpha, Iz_over_Kz)

    if t <= T:
        return (beta * (1. - (d / Z))) + (Iz_over_Kz * R_fn(t_star))
    else:
        return (beta * (1. - (d / Z))) + (Iz_over_Kz * (R_fn(t_star) - R_fn(t_star - T_star)))


def Iverson_Fig_7(t, T, Do, alpha, Iz_over_Kz):

    D = np.linspace(0.1, 6., 1000)

    x = []

    for d in D:
        x.append(Iverson_Eq_27ab(days_to_secs(6.), days_to_secs(10.), 0.000001, math.radians(15.), d))

    plt.gca().invert_yaxis()
    plt.plot(x, D)
    plt.show()


def Iverson_Fig_5(T_star):

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
    plt.xlabel('Normalised Time, $T^*$')
    plt.ylabel('R')
    legend = plt.legend()
    legend.get_frame().set_linewidth(0.)
    plt.tight_layout()
    plt.savefig('Fig_5.png')

Iverson_Fig_5(10)
Iverson_Fig_7(days_to_secs(6.), days_to_secs(10.), 0.000001, math.radians(15.), 1.)
