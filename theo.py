import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d

files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\*.txt')
files = files[:-2]
element = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\NIST\*")
X = [[], [], [], [], [], [], [], [], []]
Y = [[], [], [], [], [], [], [], [], []]
int = np.zeros(len(files))

for i in range(len(files)):
    f = open(files[i], 'r')
    for line in f:
        row = line.split()
        X[i].append(float(row[0]))
        Y[i].append(float(row[1]))
    X[i] = np.array(X[i])
    Y[i] = np.array(Y[i])

    plt.plot(X[i], Y[i])
    plt.title(files[i])
    plt.show()

Ey = [[], [], [], [], [], [], []]
Ex = [[], [], [], [], [], [], []]
F = []
for i in range(len(element)):
    f = open(element[i], 'r')
    for line in f:
        row = line.split()
        Ex[i].append(float(row[0]))
        Ey[i].append(float(row[-1]))
    Ex[i] = np.array(Ex[i])
    Ey[i] = np.array(Ey[i])

    F.append(interp1d(Ex[i], Ey[i], fill_value='extrapolate'))

    plt.plot(Ex[i], F[i](Ex[i]), 'ro')
    plt.plot(Ex[i], Ey[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.title(element[i])
    plt.grid()
    plt.show()



fCa = 20/68
fS = 16/68
fO4 = 32/86

fLi = 6/82
fB = 20/82
fO7 = 56/82

# E = np.array([50, 60, 80, 100, 150])
eff = np.array([7.01, 14.9, 30, 45, 60.1, 75, 90.1, 105, 120])
muB = F[0](eff*10**(-3))
muCa = F[1](eff*10**(-3))
muH = F[2](eff*10**(-3))
muLi = F[3](eff*10**(-3))
muO = F[4](eff*10**(-3))
muS = F[5](eff*10**(-3))
muw = F[6](eff*10**(-3))


# for i in range(6):
#     plt.plot(Ex[i], Ey[i])
#     plt.plot(eff*10**(-3), F[i](eff*10**(-3)), 'o')
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.title(element[i])
#     plt.grid()
#     plt.show()



mu_CaSO4 = muCa*fCa + muS*fS + muO*fO4
mu_LTB = muLi*fLi + muB*fB + muO*fO7
comp_w = (2/10)*muH + (8/10)*muO


a, b, c, d = np.polyfit(eff, mu_CaSO4/mu_LTB, 3)

plt.plot(eff, mu_CaSO4/muw, 'o')
plt.plot(eff, mu_LTB/muw, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.xlabel('keV')
plt.ylabel('(mu_en/ro)/(mu_en/ro)_w')
plt.title('mass energy absorption coefficients')
plt.show()
plt.plot(eff, mu_CaSO4/mu_LTB, 'o')
plt.grid()
plt.xlabel('keV')
plt.ylabel('Ratio')
plt.show()


#prøver spekter integrasjon istede for å legge sammen med vekt forholdene.

def spec_int(x, y, mu):
    D = 0
    dw = 0
    for i in range(len(y) - 1):
        D += y[i + 1]*(mu[i])*(x[i + 1] - x[i])
        dw += y[i + 1]*(F[6](x[i]))*(x[i + 1] - x[i])
    return D/dw


D_CaSO4 = np.zeros(9)
D_LTB = np.zeros(9)

# print(spec_int(X[0], Y[0], F[1](X[0])*fCa + F[4](X[0])*fS + F[3](X[0])*fO4))

for i in range(len(Y)):
    D_CaSO4[i] = spec_int(X[i]*10**(-3), Y[i], F[1](X[i]*10**(-3))*fCa + F[5](X[i]*10**(-3))*fS + F[4](X[i]*10**(-3))*fO4)
    D_LTB[i] = spec_int(X[i]*10**(-3), Y[i], F[3](X[i]*10**(-3))*fLi + F[0](X[i]*10**(-3))*fB + F[4](X[i]*10**(-3))*fO7)

    # plt.plot(Ex[0], F[1](Ex[0])*fCa + F[5](Ex[0])*fS + F[4](Ex[0])*fO4)
    # plt.plot(X[i]*10**(-3), F[1](X[i]*10**(-3))*fCa + F[5](X[i]*10**(-3))*fS + F[4](X[i]*10**(-3))*fO4, 'o')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.grid()
    # plt.title('CaSO4')
    # plt.show()
    # plt.plot(Ex[0], F[3](Ex[0])*fLi + F[0](Ex[0])*fB + F[4](Ex[0])*fO7)
    # plt.plot(X[i]*10**(-3), F[3](X[i]*10**(-3))*fLi + F[0](X[i]*10**(-3))*fB + F[4](X[i]*10**(-3))*fO7, 'o')
    # plt.yscale('log')
    # plt.xscale('log')
    # plt.grid()
    # plt.title('LTB')
    # plt.show()

plt.plot(eff, (mu_CaSO4/mu_LTB)/(np.average(mu_CaSO4/mu_LTB)), 'o')
plt.plot(eff, (D_CaSO4/D_LTB)/(np.average(D_CaSO4/D_LTB)), 'o')
plt.legend(['(mu_en/ro) ratio', 'Spectral average method'])
plt.grid()
plt.xlabel('keV')
plt.ylabel('Normalized ratio')
# plt.title('Ratio ')
plt.show()

plt.plot(eff, D_CaSO4, 'o')
plt.plot(eff, D_LTB, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.xlabel('keV')
plt.ylabel('mu_en/ro)/(mu_en/ro)_w')
plt.title('Spectral average')
plt.show()
plt.plot(eff, D_CaSO4/ D_LTB, 'o')
plt.grid()
plt.xlabel('keV')
plt.ylabel('ratio')
plt.title('Spectral average')
plt.show()



HVL_Al = [0.243, 1.22, 2.11, 5.83, 10.7, 12.7, 15.0, 16.6] #mm
HVL_Cu = [0.00814, 0.0378, 0.0698, 0.264, 0.730, 1.17, 1.79, 2.39]

# plt.plot(HVL_Al, D_CaSO4/D_LTB, 'o')
# plt.grid()
# plt.xlabel('HVL_Al [mm]')
# plt.ylabel('Ratio')
# plt.title('Ratio from x - ray spectra')
# plt.show()
# plt.plot(HVL_Cu, D_CaSO4/D_LTB, 'o')
# plt.grid()
# plt.xlabel('HVL_Cu [mm]')
# plt.ylabel('Ratio')
# plt.title('Ratio from x - ray spectra')
# plt.show()


plt.plot(Ex[0], F[1](Ex[0])*fCa + F[5](Ex[0])*fS + F[4](Ex[0])*fO4)
plt.plot(Ex[0], F[3](Ex[0])*fLi + F[0](Ex[0])*fB + F[4](Ex[0])*fO7)
plt.plot(Ex[0], F[2](Ex[0])*(2/10) + F[4](Ex[0])*(8/10))
plt.legend(['CaSO4', 'LTB', 'H2O'])
plt.xscale('log')
plt.yscale('log')
plt.title('composite mu_en/ro')
plt.grid()
plt.xlabel('MeV')
plt.ylabel('mu_en/ro [cm^2/g]')
plt.show()

plt.plot(eff, comp_w)
plt.plot(eff, muw)
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlabel('MeV')
plt.ylabel('mu_en/ro [cm^2/g]')
plt.show()

# print('Co60 CaSO4 og LTB')
# print(F[1](1.068253259235031)*fCa + F[5](1.068253259235031)*fS + F[4](1.068253259235031)*fO4, F[3](1.068253259235031)*fLi + F[0](1.068253259235031)*fB + F[4](1.068253259235031)*fO7)
# print('Blindstudie oppsett CaSO4 og LTB')
# print(F[1](0.0974)*fCa + F[5](0.0974)*fS + F[4](0.0974)*fO4, F[3](0.0974)*fLi + F[0](0.0974)*fB + F[4](0.0974)*fO7)
