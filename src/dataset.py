from pathlib import Path
import numpy as np

DATADIR = Path("data")
assert DATADIR.exists()


PROTO_FILENAME_LES = "WT{:03d}.T1200s.ufluct"

alpha_zero = 1.69292
R = 46.5
center = (463.136, 64)
Uamb = 8.0  # m/s

n1 = 10325
n2 = 68
n3 = 68

v_shear = [
    0.07273814,
    0.20838735,
    0.33097443,
    0.44049183,
    0.5369399,
    0.6203186,
    0.6906278,
    0.74786747,
    0.79203767,
    0.82313824,
    0.8411693,
    0.8519544,
    0.861958,
    0.8712956,
    0.8800561,
    0.8883115,
    0.8961209,
    0.9035334,
    0.9105902,
    0.91732615,
    0.9237713,
    0.92995155,
    0.93588936,
    0.9416044,
    0.94711393,
    0.95243335,
    0.95757633,
    0.96255505,
    0.9673804,
    0.97206223,
    0.9766096,
    0.98103046,
    0.9853322,
    0.9895216,
    0.99360484,
    0.99758744,
    1.0014747,
    1.0052716,
    1.0089823,
    1.0126109,
    1.0161614,
    1.0196375,
    1.0230421,
    1.0263784,
    1.0296495,
    1.032858,
    1.0360065,
    1.0390971,
    1.0421324,
    1.0451143,
    1.0480448,
    1.0509257,
    1.0537591,
    1.0565463,
    1.0592892,
    1.0619891,
    1.0646474,
    1.0672656,
    1.069845,
    1.0723867,
    1.0748919,
    1.0773618,
    1.0797973,
    1.0821996,
    1.0845696,
    1.0869081,
    1.0892161,
    1.0914946,
    1.0937442,
    1.0959656,
    1.0981598,
    1.1003274,
    1.102469,
    1.1045853,
    1.106677,
    1.1087449,
    1.1107892,
    1.1128107,
    1.1148099,
    1.1167873,
    1.1187434,
    1.1206787,
    1.1225938,
    1.1244888,
    1.1263646,
    1.1282214,
    1.1300595,
    1.1318794,
    1.1336817,
    1.1354663,
    1.1372341,
    1.138985,
    1.1407195,
    1.142438,
    1.1441408,
    1.1458282,
    1.1475005,
    1.1491579,
    1.1508007,
    1.1524293,
    1.1540438,
    1.1556447,
    1.1572319,
    1.158806,
    1.1603669,
    1.1619151,
    1.1634507,
    1.164974,
    1.1664852,
    1.1679844,
    1.1694719,
    1.1709478,
    1.1724124,
    1.1738659,
    1.1753083,
    1.17674,
    1.1781611,
    1.1795717,
    1.1809721,
    1.1823623,
    1.1837425,
    1.185113,
    1.1864737,
    1.187825,
    1.1891668,
    1.1904994,
    1.1918229,
    1.1931375,
    1.1944432,
    1.1957402,
    1.1970286,
    1.1983086,
    1.1995802,
    1.2008436,
    1.2020988,
    1.2033461,
    1.2045856,
    1.2058171,
    1.207041,
    1.2082573,
    1.2094662,
    1.2106677,
    1.211862,
    1.2130489,
    1.2142289,
    1.2154018,
    1.2165679,
    1.2177271,
    1.2188795,
    1.2200253,
    1.2211645,
    1.2222972,
    1.2234236,
    1.2245436,
    1.2256573,
    1.2267648,
    1.2278663,
    1.2289617,
    1.2300512,
    1.2311347,
    1.2322124,
    1.2332844,
    1.2343506,
    1.2354113,
    1.2364663,
    1.2375159,
    1.2385601,
    1.2395988,
    1.2406322,
    1.2416604,
    1.2426833,
    1.2438128,
    1.2451822,
    1.2468402,
    1.2488449,
    1.251264,
    1.2541767,
    1.2576747,
    1.2618628,
    1.2668589,
    1.2727945,
    1.2798127,
    1.2880661,
    1.2977142,
    1.3089174,
    1.3218323,
    1.3366059,
    1.3533686,
    1.3722295,
    1.3932723,
    1.4165531,
    1.4421021,
]


def get_LES_coords():
    """
    Lateral and vertical coordinates of LES field as described by data set
    documentation.
    """
    n1 = 10325
    n2 = 68
    n3 = 68

    L1 = 9601.319092  # [m]
    L2 = 118.310131  # [m]
    L3 = 118.425435  # [m]
    x = np.linspace(0, L2, n2) - L2 / 2
    y = np.linspace(0, L3, n3) - L3 / 2

    return x, y


def get_LES_ref():

    fn = DATADIR / "FreeStream.T1200s.ufluct"
    LES = np.fromfile(fn, dtype=np.single)
    LES = np.reshape(LES, (n1, n2, n3))
    LES = np.swapaxes(LES, 1, 2)
    LES = LES * 8

    return LES


def get_LES_field(R):
    """
    Load a wind field as a 3D numpy array  with coordinates: (longitudinal, lateral, vertical).
    """
    fn = DATADIR / PROTO_FILENAME_LES.format(int(R))
    LES = np.fromfile(fn, dtype=np.single)
    LES = np.reshape(LES, (n1, n2, n3))
    LES = np.swapaxes(LES, 1, 2)
    # LES = np.flipud(LES)
    # LES -= np.tile(v_shear, (192, 1)).T
    LES = LES * 8

    return LES


if __name__ == "__main__":
    field = get_LES_field(1)
    ref = get_LES_ref()

    field = field - ref
    import matplotlib.pyplot as plt

    plt.imshow(ref[1000, :, :])
    plt.savefig("asdf.png")

    asdf = field.mean(axis=0).mean(axis=1)
    plt.close()
    plt.plot(asdf)
    plt.savefig("asdf2.png")
