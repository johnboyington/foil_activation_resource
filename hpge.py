import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


# each efficiency calibration will be stored within this dictionary
hpge_efficiency = {}


def ksu(erg):
    '''
    This uses data collected from a ksu hpge efficiency calculation and fits a curve to it
    '''
    def efficiency(erg, a, b, c, d, e):
        erg = np.log(erg)
        eff = a + b * erg + c * erg**2 + d * erg**3 + e * erg**4
        return np.exp(eff)

    x = np.array([60, 86.5, 105.3, 176.3, 427.90, 600.6, 834.8, 1173.2, 1332.5])
    y = np.array([2.39e-2, 7.5e-2, 9.27e-2, 9.35e-2, 4.77e-2, 3.51e-2, 2.64e-2, 1.87e-2, 1.69e-2])
    p0 = [-143.191, 91.147, -21.757, 2.279, -0.089]
    coeffs = curve_fit(efficiency, x, y, p0)
    coeffs = -2.48945787e+02, 1.67284939e+02, -4.21090055e+01, 4.67244596e+00, -1.93763734e-01
    return efficiency(erg, *coeffs)


hpge_efficiency['ksu'] = ksu


def perfect(erg):
    '''
    This is a theoretically perfect detector
    '''
    return 1.0


hpge_efficiency['perfect'] = perfect


def wisconsin(E):
    C1 = -4.314E+01
    C2 = 3.835E+01
    C3 = -7.515E+00
    C4 = 4.362E-01
    A = 8.03E-03
    DI = 5.310  # cm
    DL = 700.000  # um
    AI = 0.000  # deg
    T1 = 500.000  # um
    det_model = A*np.exp((-T1*U1-T2*U2-T3*U3-T4*U4-DL*U5) / np.cos(AI))*(1-np.exp(-DI*U5/np.cos(AI)))
    return det_model * (C1 + C2*np.log(E) + C3*np.log(E)^2 + C4*np.log(E)^3)


hpge_efficiency['wisconsin'] = wisconsin










if __name__ == '__main__':
    x = np.array([60, 86.5, 105.3, 176.3, 427.90, 600.6, 834.8, 1173.2, 1332.5])
    y = np.array([2.39e-2, 7.5e-2, 9.27e-2, 9.35e-2, 4.77e-2, 3.51e-2, 2.64e-2, 1.87e-2, 1.69e-2])
    x1 = np.linspace(0, 2000, 2000)
    y1 = hpge_efficiency['ksu'](x1)
    plt.plot(x1, y1, 'b', x, y, 'ko')
    plt.ylim(0)
