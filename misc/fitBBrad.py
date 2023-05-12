import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# bbcon = 14387768.8  # nm K
steps = 500


def g(x, a, mu, s):
    return a * np.exp(-0.5*((x-mu)/s)**2)


def gn(x, *args):
    # assume len(args) is divisible by 3
    n = len(args) // 3
    y = x * 0
    for i in range(n):
        y += g(x, *args[3*i:3*i+3])
    return y


def bb(xT):
    return 1 / (xT**5 * np.expm1(1/xT))


x = np.linspace(1/steps, 3, steps)
y = bb(x)

# plt.figure(0)
# plt.plot(x, y, label="bb(x)")
popt0 = [20, 5, 5]
popt = []
plt.figure()
for i in range(6):
    if i == 0:
        print(f"Fit 1 mode:")
    else:
        print(f"Fit {i+1} modes:")
    popt, pcov = curve_fit(gn, x, y, p0=list(popt) + popt0, maxfev=1000000)
    for j in range(i+1):
        a = popt[3*j]
        m = popt[3*j+1]
        s = popt[3*j+2]
        print(f"max_a{j+1}: {a * np.sqrt(2*np.pi) * s},",
              f"mu{j+1}: {m}, s{j+1}: {s}")
    yn = gn(x, *popt)
    d2 = (yn-y)**2
    m_d = np.sqrt(max(d2))
    t_d2 = np.sum(d2)
    print(f"max_diff: {m_d}, tot_d2: {t_d2}")
    plt.subplot(2, 3, i+1)
    plt.plot(x, y, label="bb(x)")
    plt.plot(x, yn, label=f"g{i+1} fit")
    plt.xlim([0, 1.5])
    plt.ylim([0, 25])
    plt.legend()

plt.show()
