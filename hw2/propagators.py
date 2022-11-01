from scipy.special import jn, yn, kn
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pylab as pylab
plt.style.use("science")  # requires SciencePlots

params = {
    'legend.fontsize': 15,
    'figure.figsize': (5, 5),
    'axes.labelsize': 15,
    'axes.titlesize': 20,
    'xtick.labelsize': 15,
    'ytick.labelsize': 15,
    'figure.titlesize': 30,
    'xtick.major.size': 4,
    'xtick.minor.size': 2,
    'xtick.major.width': 1,
    'xtick.minor.width': 1,
    'ytick.major.size': 4,
    'ytick.minor.size': 2,
    'ytick.major.width': 1,
    'ytick.minor.width': 1,
    'font.size': 20 # for annotate
}
pylab.rcParams.update(params)

N = 1000
m = 1
T = 10
A = 0.05
x = np.linspace(0.1, T, N)
pi = np.pi

def phi_spacelike(x, m=m):
    return m / 4 / pi**2 / x * kn(1, x * m)

def phi_timelike(x, m=m):
    return - 1j * m / 4 / pi / x * (jn(1, m * x) - 1j * yn(1, m * x))

def pi_spacelike(x, m=m):
    return - m**2 / 4 / pi**2 / x**2 * kn(2, m * x)

def pi_timelike(x, m=m):
    return 1j * m**3 / 4 / pi / x * (jn(1, m * x) - 1j * yn(1, m * x) - 1/m/x * jn(2, m * x) + 1j / m / x * yn(2, m*x))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
ax1.set_title(r"$D(x - y)$")
ax1.plot(x, phi_spacelike(x), color="r", label=r"$(x - y)^2 < 0$")
ax1.plot(x, np.real(phi_timelike(x)), color="b", label=r"$(x - y)^2 > 0$ ($\Re$)")
ax1.plot(x, np.imag(phi_timelike(x)), ls="--", color="b", label=r" $(x - y)^2 > 0$ ($\Im$)")
ax2.set_title(r"$D_{\pi}(x - y)$")
ax2.plot(x, pi_spacelike(x), ls="-", color="r")
ax2.plot(x, np.real(pi_timelike(x)), ls="-", color="b")
ax2.plot(x, np.imag(pi_timelike(x)), ls="--", color="b")
ax1.set_ylim(-A, A)
ax1.set_xlim(0, T)
ax2.set_xlim(0, T)

xticks = [1/m]
labels = [r"$\dfrac{1}{m}$"]
ax1.set_xticks(xticks)
ax2.set_xticks(xticks)
ax1.set_xticklabels(labels)
ax2.set_xticklabels(labels)

yticks = [0, A]
ax1.set_yticks(yticks)
labels = ["0", f"{A:.2f}"] 
ax1.set_yticklabels(labels)

ax1.legend()
ax1.axhline(0, color="k", linestyle="--")
ax2.axhline(0, color="k", linestyle="--")
plt.subplots_adjust(hspace=0, wspace=0)
plt.savefig("amplitudes.pdf", bbox_inches="tight")
plt.show()


# print(phi_spacelike(1, 1))
# print(phi_spacelike(2, 1))
# print(phi_spacelike(np.e, 1))

# print(phi_spacelike(1/2, 2))
# print(phi_spacelike(1, 2))
# print(phi_spacelike(3/2, 2))
    
