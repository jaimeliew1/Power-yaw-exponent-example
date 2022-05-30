from pathlib import Path

import numpy as np
import pandas as pd
from scipy import interpolate, optimize

import dataset
from lib import ellipse_coords, yaw_transform

OUTDIR = Path("out")
OUTDIR.mkdir(parents=True, exist_ok=True)

rotor_D = 2 * dataset.R
Yaws = [-30, -20, -10, 0, 10, 20, 30]
Rs = np.arange(2, 24, 1)


def get_radial_windspeed(Rs, X, Y, field, center=(0.0, 0.0)):
    """
    Eq. 12. Converts 2D wind field in cartesian coordinates to radial
    coordinates with azimuthal averaging.
    """
    field_f = interpolate.interp2d(X, Y, field)

    U = np.zeros(len(Rs))
    for i, R in enumerate(Rs):
        circx, circy, theta = ellipse_coords(center[0], center[1], R)
        circx, circy = (
            circx[abs(theta - np.pi) > np.deg2rad(30)],
            circy[abs(theta - np.pi) > np.deg2rad(30)],
        )
        U[i] = np.mean([field_f(x, y)[0] for x, y in zip(circx, circy)])

    return U


def power(r, U, yaw, alpha_zero):
    """
    Eq. 8. Returns the power output of a yawed turbine with yaw-transformed wind
    speed profile in radial coordinates.
    """
    r, U = r[r <= rotor_D / 2], U[r <= rotor_D / 2]
    P = np.cos(np.deg2rad(yaw)) ** alpha_zero * np.trapz(r * U**3, r)

    return P


def _fit_func(x, alpha):
    """x in degrees"""
    return np.cos(np.deg2rad(x)) ** alpha


def fit_alpha(yaw, P_norm):
    """
    Minimization described in Eq. 10. Returns the value of alpha for a set of
    non_dimensional_power-yaw value pairs.
    """
    popt, _ = optimize.curve_fit(_fit_func, yaw, P_norm)
    return popt[0]


if __name__ == "__main__":
    Alphas = np.zeros(len(Rs))
    yawDP_table = []  # placeholder for yaw-distance-power table
    X, Y = dataset.get_LES_coords()
    for k, R in enumerate(Rs):
        field = dataset.get_LES_field(R) - dataset.get_LES_ref()

        # Average over time
        field = field.mean(axis=0)

        # get average wind speed over each annulus, and interpolate
        Rs = np.linspace(0, rotor_D, 100)
        Udef = get_radial_windspeed(Rs, X, Y, field, center=dataset.center)
        Udef_func = interpolate.UnivariateSpline(Rs, Udef, s=0)

        # Get reference power output at yaw = 0.
        P0 = power(
            Rs, yaw_transform(Rs, Udef_func, 0) + dataset.Uamb, 0, dataset.alpha_zero
        )

        # Transform wake deficit, and calculate relative power
        P = np.zeros(len(Yaws))
        for j, yaw in enumerate(Yaws):
            Udef_trans = yaw_transform(Rs, Udef_func, yaw)
            P[j] = power(Rs, Udef_trans + dataset.Uamb, yaw, dataset.alpha_zero) / P0
            yawDP_table.append([yaw, R, P[j]])

        # Fit alpha
        Alphas[k] = fit_alpha(Yaws, P)
        print("R = {}, alpha={}".format(R, Alphas[k]))
        print(P)

    yawDP_table = pd.DataFrame(yawDP_table, columns=["yaw", "R", "P"])
    yawDP_table.to_csv(OUTDIR / "yawDP_analytical_LES.csv")

    # write alphas to file
    with open(OUTDIR / "alpha_analytical_LES.txt", "w") as f:
        for d, a in zip(Rs, Alphas):
            f.write(f"{d}, {a}\n")
