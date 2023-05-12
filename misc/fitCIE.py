import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
import pandas

data = pandas.read_csv("CIE_xyz_1931_2deg.csv",
                       names=['lambda', 'X', 'Y', 'Z'])

lamda = data['lambda'].to_numpy(dtype=np.float32)
X = data['X'].to_numpy(dtype=np.float32)
Y = data['Y'].to_numpy(dtype=np.float32)
Z = data['Z'].to_numpy(dtype=np.float32)

steps = len(lamda)


def g(x, a, mu, s):
    return a * np.exp(-0.5*((x-mu)/s)**2)


def gn(x, *args):
    # assume len(args) is divisible by 3
    n = len(args) // 3
    y = np.zeros_like(x, dtype=np.float32)
    for i in range(n):
        y += g(x, *args[3*i:3*i+3])
    return y


# plt.figure(0)
# plt.plot(lamda, X, label="X")
# plt.plot(lamda, Y, label="Y")
# plt.plot(lamda, Z, label="Z")
# plt.legend()
# plt.show()

poptX = [1.07, 601, 80,
         0.345, 448, 30,
         0.1, 550, 50]
poptY = [1, 554, 60,
         0.2, 600, 50,
         0.2, 500, 50]
poptZ = [1.778, 445, 40,
         0.1, 500, 30,
         0.1, 450, 100]

plt.figure()
for idx, (y, name, popt0) in enumerate([(X, "X", poptX),
                                       (Y, "Y", poptY),
                                       (Z, "Z", poptZ)]):
    for i in range(3):
        print(f"Fitting {i+1} mode(s) to CIE {name}:")
        popt, pcov = curve_fit(
            gn, lamda, y, p0=popt0[:3*(i+1)], method="dogbox", maxfev=100000)
        tot_area = 0
        for j in range(i+1):
            a = popt[3*j]
            m = popt[3*j+1]
            s = popt[3*j+2]
            area = a * np.sqrt(2*np.pi) * s
            tot_area += area
            print(f"area{j+1}: {area}, a{j+1}: {a}",
                  f"mu{j+1}: {m}, s{j+1}: {s}")
        print(f"tot_area: {tot_area}")

        yn = gn(lamda, *popt)
        d2 = (yn-y)**2
        m_d = np.sqrt(max(d2))
        t_d2 = np.sum(d2)
        print(f"max_diff: {m_d}, tot_d2: {t_d2}")
        print("---------")

        # plt.figure()
        plt.subplot2grid((3, 3), (idx, i))
        plt.plot(lamda, y, label=name)
        plt.plot(lamda, yn, label=f"{name} fit, {i+1} mode(s)")
        plt.legend()

plt.show()
