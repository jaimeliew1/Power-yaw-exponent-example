"""
A collection of useful functions for this paper.
"""
import numpy as np
from scipy import interpolate


def ellipse_r(theta, R, e=0):
    # return R**2*e/np.sqrt(R**2*e**2*cos(theta)**2 + R**2*sin(theta)**2)
    return R * np.sqrt(1 - e**2) / np.sqrt(1 - e**2 * np.cos(theta) ** 2)


def CDF(r, R, yaw):
    """
    Calculates the cumulative density function of an elliptic path as described
    in Eq. 5 in Liew, J (2020).

    arguments:
        r (float or numpy array): radial position in meteorological coordinates [m]
        R (float): radial position in rotor coordinates [m]
        yaw (float): yaw angle of rotor [deg]
    """
    yaw = np.deg2rad(yaw)
    return 1 - 2 / np.pi * np.arccos(
        np.sqrt((r**2 - R**2 * np.cos(yaw) ** 2) / (r**2 * np.sin(yaw) ** 2))
    )


def PDF(r, R, yaw):
    """
    Calculates the probability density function of an elliptic path as described
    in Eq. 5 in Liew, J (2020).

    arguments:
        r (float or numpy array): radial position in meteorological coordinates [m]
        R (float): radial position in rotor coordinates [m]
        yaw (float): yaw angle of rotor [deg]
    """
    e = abs(np.sin(np.deg2rad(yaw)))
    out = (
        2
        / np.pi
        * (1 - e**2)
        * R**2
        / (
            r
            * np.sqrt(
                (1 - e**2) * (R**2 - r**2) * (r**2 - (1 - e**2) * R**2)
            )
        )
    )
    out[r <= R * np.sqrt(1 - e**2)] = 0
    out[r >= R] = 0
    return out


def rad2cart_deficit(R, Udef, x0=0, y0=0):
    interp_func = interpolate.interp1d(R, Udef, bounds_error=False)
    return lambda x, y: interp_func(np.sqrt((x - x0) ** 2 + (y - y0) ** 2))


def ellipse_coords(x0, y0, R, e=0):
    # returns x, y, theta coords for an ellipse. todo: proper definition of eccentricity
    theta = np.linspace(0, 2 * np.pi, 500)
    circx, circy = R * np.sqrt(1 - e**2) * np.sin(theta) + x0, R * np.cos(theta) + y0
    return circx, circy, theta


def yaw_transform(Rs, Ufunc, yaw, eps=1e-3, disc=100):
    """
    Transforms an assymetric radial windspeed function in meteorological
    coordinates into rotor coordinates as described in Eq. 4 in Liew, J (2020).

    Arguments
    Rs: list of floats
        Radial positions along blade to consider.
    Ufunc: scipy.interpolate.UnivariateSpline
        Radial wind speed function. Takes radial position as argument.
    yaw: float
        Yaw angle of rotor [deg]
    eps: float (default: 1e-3)
        A small number added/subracted to Rmin and Rmax to prevent arithmetic error.
    disc: int (default: 100)
        Number of discretisations used to perform trapezoidal integration.

    Returns
    -------
    out: list of floats
        Transformed wind speeds at radial locations defined in Rs.
    """
    if yaw == 0:
        return Ufunc(Rs)
    out = np.zeros(len(Rs))

    for i, r in enumerate(Rs):
        # Bounds of CDF in Eq. 5
        Rmin, Rmax = r * np.cos(np.deg2rad(yaw)), r
        rs = np.linspace(Rmin + eps, Rmax - eps, disc)

        # derivative of U with respect to radius. Used in Eq. 4
        dudr = Ufunc.derivative()(rs)
        cdf = CDF(rs, r, yaw)

        # Eq. 4.
        out[i] = Ufunc(Rmax) - np.trapz(dudr * cdf, rs)
    out[np.isnan(out)] = 0

    return out
