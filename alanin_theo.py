import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d

files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Alanin spek\*.txt')
element = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\NIST 2\*")
X = [[], [], [], [], []]
Y = [[], [], [], [], []]
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

Ey = [[], [], []]
Ex = [[], [], []]
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

    plt.plot(Ex[0], F[i](Ex[0]), 'ro')
    plt.plot(Ex[i], Ey[i])
    plt.yscale('log')
    plt.xscale('log')
    plt.title(element[i])
    plt.grid()
    plt.show()


plt.plot(Ex[0], F[0](Ex[0]))
plt.plot(Ex[0], F[1](Ex[0]))
plt.plot(Ex[0], F[2](Ex[0]))
plt.title('Mass absorption coef.')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend(['Alanine', 'Cu', 'Water'])
plt.xlabel('MeV')
plt.ylabel('(mu_en/rho)')
plt.show()

eff = np.linspace(0.01, 0.2, 100)



# plt.plot(eff*10**3, F[2](eff)/(F[2](eff)))
# plt.plot(eff*10**3, F[0](eff)/(F[2](eff)))
# plt.plot(eff*10**3, (0.995*F[0](eff) + 0.005*F[1](eff))/(F[2](eff)))
# plt.plot(eff*10**3, (0.99*F[0](eff) + 0.01*F[1](eff))/(F[2](eff)))
# plt.plot(eff*10**3, (0.985*F[0](eff) + 0.015*F[1](eff))/(F[2](eff)))
# plt.plot(eff*10**3, (0.98*F[0](eff) + 0.02*F[1](eff))/(F[2](eff)))
# plt.plot(eff*10**3, (0.975*F[0](eff) + 0.025*F[1](eff))/(F[2](eff)))
# plt.plot(eff*10**3, (0.97*F[0](eff) + 0.03*F[1](eff))/(F[2](eff)))
# plt.xscale('log')
# plt.grid()
# plt.legend(['Water', 'Pure Alanine', '0.5% Cu', '1.0% Cu', '1.5% Cu', '2.0% Cu', '2.5% Cu', '3.0% Cu'])
# plt.xlabel('keV')
# plt.show()


#Regner ut med prosenter basert på Ravis beregning av effektiv atomvekt

plt.plot(eff*10**3, F[2](eff)/(F[2](eff)))
plt.plot(eff*10**3, F[0](eff)/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.00355)*F[0](eff) + 0.00355*F[1](eff))/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.007081976)*F[0](eff) + 0.007081976*F[1](eff))/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.01058548)*F[0](eff) + 0.01058548*F[1](eff))/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.014064348)*F[0](eff) + 0.014064348*F[1](eff))/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.017518837)*F[0](eff) + 0.017518837*F[1](eff))/(F[2](eff)))
plt.plot(eff*10**3, ((1 - 0.020949203)*F[0](eff) + 0.020949203*F[1](eff))/(F[2](eff)))
# plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.legend(['Water', 'Pure Alanine', '0.5 mol% Cu', '1.0 mol% Cu', '1.5 mol% Cu', '2.0 mol% Cu', '2.5 mol% Cu', '3.0 mol% Cu'])
plt.xlabel('keV')
plt.ylabel('(mu_en/rho)Cu, Alanin/ (mu_en/rho)Water')
plt.show()



test = np.linspace(0.025, 0.120, 5)
plt.plot(test*10**3, F[2](test)/(F[2](test)))
plt.plot(test*10**3, F[0](test)/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.00355)*F[0](test) + 0.00355*F[1](test))/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.007081976)*F[0](test) + 0.007081976*F[1](test))/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.01058548)*F[0](test) + 0.01058548*F[1](test))/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.014064348)*F[0](test) + 0.014064348*F[1](test))/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.017518837)*F[0](test) + 0.017518837*F[1](test))/(F[2](test)), 'o')
plt.plot(test*10**3, ((1 - 0.020949203)*F[0](test) + 0.020949203*F[1](test))/(F[2](test)), 'o')
# plt.xscale('log')
# plt.yscale('log')
plt.grid()
plt.legend(['Water', 'Pure Alanine', '0.5 mol% Cu', '1.0 mol% Cu', '1.5 mol% Cu', '2.0 mol% Cu', '2.5 mol% Cu', '3.0 mol% Cu'])
plt.xlabel('keV')
plt.ylabel('(mu_en/rho)Cu, Alanin/ (mu_en/rho)Water')
plt.show()



# #prøver spekter integrasjon istede for å legge sammen med vekt forholdene.

def spec_int(x, y, mu):
    D = 0
    dw = 0
    for i in range(len(y) - 1):
        D += y[i + 1]*(mu[i])*(x[i + 1] - x[i])
        dw += y[i + 1]*(F[2](x[i]))*(x[i + 1] - x[i])
    return D/dw



pure = np.zeros(5)
Cu05 = np.zeros(5)
Cu1 = np.zeros(5)
Cu15 = np.zeros(5)
Cu2 = np.zeros(5)
Cu25 = np.zeros(5)
Cu3 = np.zeros(5)


for i in range(len(Y)):
    pure[i] = spec_int(X[i]*10**(-3), Y[i], F[0](X[i]*10**(-3)))
    Cu05[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.00355)*F[0](X[i]*10**(-3)) + 0.00355*F[1](X[i]*10**(-3))))
    Cu1[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.007081976)*F[0](X[i]*10**(-3)) + 0.007081976*F[1](X[i]*10**(-3))))
    Cu15[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.01058548)*F[0](X[i]*10**(-3)) + 0.01058548*F[1](X[i]*10**(-3))))
    Cu2[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.014064348)*F[0](X[i]*10**(-3)) + 0.014064348*F[1](X[i]*10**(-3))))
    Cu25[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.017518837)*F[0](X[i]*10**(-3)) + 0.017518837*F[1](X[i]*10**(-3))))
    Cu3[i] = spec_int(X[i]*10**(-3), Y[i], ((1 - 0.020949203)*F[0](X[i]*10**(-3)) + 0.020949203*F[1](X[i]*10**(-3))))





plt.plot(test*10**3, F[2](test)/(F[2](test)))
plt.plot(test*10**3, pure, '-o')
plt.plot(test*10**3, Cu05, '-o')
plt.plot(test*10**3, Cu1, '-o')
plt.plot(test*10**3, Cu15, '-o')
plt.plot(test*10**3, Cu2, '-o')
plt.plot(test*10**3, Cu25, '-o')
plt.plot(test*10**3, Cu3, '-o')



# plt.plot(test*10**3, F[0](test)/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.00355)*F[0](test) + 0.00355*F[1](test))/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.007081976)*F[0](test) + 0.007081976*F[1](test))/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.01058548)*F[0](test) + 0.01058548*F[1](test))/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.014064348)*F[0](test) + 0.014064348*F[1](test))/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.017518837)*F[0](test) + 0.017518837*F[1](test))/(F[2](test)))
# plt.plot(test*10**3, ((1 - 0.020949203)*F[0](test) + 0.020949203*F[1](test))/(F[2](test)))
plt.title('Spectral average')
plt.grid()
plt.legend(['Water', 'Pure Alanine', '0.5 mol% Cu', '1.0 mol% Cu', '1.5 mol% Cu', '2.0 mol% Cu', '2.5 mol% Cu', '3.0 mol% Cu'], bbox_to_anchor=(1.0, 1.0))
plt.xlabel('keV')
plt.ylabel('(mu_en/rho)Cu, Alanin/ (mu_en/rho)Water')
plt.show()









#     # plt.plot(Ex[0], F[1](Ex[0])*fCa + F[5](Ex[0])*fS + F[4](Ex[0])*fO4)
#     # plt.plot(X[i]*10**(-3), F[1](X[i]*10**(-3))*fCa + F[5](X[i]*10**(-3))*fS + F[4](X[i]*10**(-3))*fO4, 'o')
#     # plt.yscale('log')
#     # plt.xscale('log')
#     # plt.grid()
#     # plt.title('CaSO4')
#     # plt.show()
#     # plt.plot(Ex[0], F[3](Ex[0])*fLi + F[0](Ex[0])*fB + F[4](Ex[0])*fO7)
#     # plt.plot(X[i]*10**(-3), F[3](X[i]*10**(-3))*fLi + F[0](X[i]*10**(-3))*fB + F[4](X[i]*10**(-3))*fO7, 'o')
#     # plt.yscale('log')
#     # plt.xscale('log')
#     # plt.grid()
#     # plt.title('LTB')
#     # plt.show()
#
# plt.plot(eff, mu_CaSO4/mu_LTB, 'o')
# plt.plot(eff, D_CaSO4/D_LTB, 'o')
# plt.legend(['(mu_en/ro) ratio', 'Spectral average'])
# plt.grid()
# plt.xlabel('keV')
# plt.ylabel('Ratio')
# plt.title('Ratio from x - ray spectra')
# plt.show()
#
# plt.plot(eff, D_CaSO4, 'o')
# plt.plot(eff, D_LTB, 'o')
# plt.legend(['CaSO4', 'LTB'])
# plt.grid()
# plt.title('Spektral integrasjon')
# plt.xlabel('keV')
# # plt.title('Ratio from x - ray spectra')
# plt.show()
#
#
# HVL_Al = [0.243, 1.22, 2.11, 5.83, 10.7, 12.7, 15.0, 16.6] #mm
# HVL_Cu = [0.00814, 0.0378, 0.0698, 0.264, 0.730, 1.17, 1.79, 2.39]
#
# # plt.plot(HVL_Al, D_CaSO4/D_LTB, 'o')
# # plt.grid()
# # plt.xlabel('HVL_Al [mm]')
# # plt.ylabel('Ratio')
# # plt.title('Ratio from x - ray spectra')
# # plt.show()
# # plt.plot(HVL_Cu, D_CaSO4/D_LTB, 'o')
# # plt.grid()
# # plt.xlabel('HVL_Cu [mm]')
# # plt.ylabel('Ratio')
# # plt.title('Ratio from x - ray spectra')
# # plt.show()
#
#
# plt.plot(Ex[0], F[1](Ex[0])*fCa + F[5](Ex[0])*fS + F[4](Ex[0])*fO4)
# plt.plot(Ex[0], F[3](Ex[0])*fLi + F[0](Ex[0])*fB + F[4](Ex[0])*fO7)
# plt.plot(Ex[0], F[2](Ex[0])*(2/10) + F[4](Ex[0])*(8/10))
# plt.legend(['CaSO4', 'LTB', 'H2O'])
# plt.xscale('log')
# plt.yscale('log')
# plt.title('composite mu_en/ro')
# plt.grid()
# plt.xlabel('MeV')
# plt.ylabel('mu_en/ro [cm^2/g]')
# plt.show()
#
# plt.plot(eff, comp_w)
# plt.plot(eff, muw)
# plt.xscale('log')
# plt.yscale('log')
# plt.grid()
# plt.xlabel('MeV')
# plt.ylabel('mu_en/ro [cm^2/g]')
# plt.show()
#
# print('Co60 CaSO4 og LTB')
# print(F[1](1.068253259235031)*fCa + F[5](1.068253259235031)*fS + F[4](1.068253259235031)*fO4, F[3](1.068253259235031)*fLi + F[0](1.068253259235031)*fB + F[4](1.068253259235031)*fO7)
# print('Blindstudie oppsett CaSO4 og LTB')
# print(F[1](0.0974)*fCa + F[5](0.0974)*fS + F[4](0.0974)*fO4, F[3](0.0974)*fLi + F[0](0.0974)*fB + F[4](0.0974)*fO7)
