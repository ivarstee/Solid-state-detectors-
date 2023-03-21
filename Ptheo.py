import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import interp1d

def Cstop(T, Z, z, A, I): #z er ladning til stråling (+1 for proton), Z er atomnummer til stoff, A er antall nukleoner til stoff, I er mean exitation potential
    m = 1.67262192*10**(-27)
    c = 2.99792*10**8
    con = 6.24150907*10**(12) # konverterer til MeV så pass på at innput energi er i denne enheten
    beta = np.sqrt(1 - (1/(1 + T/(con*(m*c**2))))**2) #beta ser riktig ut
    # print((1/(np.sqrt(1 - beta**2)) - 1)*con*(m*c**2))
    # print(13.8373 + np.log(beta**2/(1 - beta**2)) - beta**2 )
    # print(0.3071*Z*z**2/(A*beta**2))
    return 0.3071*Z*z**2/(A*beta**2)*(13.8373 + np.log(beta**2/(1 - beta**2)) - beta**2 - np.log(I))


#finner data om stoffene i Attix tabell B1 og B2, I finner jeg fra NIST: https://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html

fCa = 20/68
ICa = 191 #eV
ACa = 40.08

fS = 16/68
IS =180
AS = 32.06

fO4 = 32/86
IO = 95.0
AO = 16.0

fLi = 6/82
ILi = 40.0
ALi = 6.941

fB = 20/82
IB = 76.0
AB = 10.81

fO7 = 56/82


files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\stopping\*.txt')
X = [[], [], [], [], [], [], [], []]
Y = [[], [], [], [], [], [], [], []]
int = np.zeros(len(files))

for i in range(len(files)):
    f = open(files[i], 'r')
    for line in f:
        row = line.split()
        X[i].append(float(row[0]))
        Y[i].append(float(row[1]))
    X[i] = np.array(X[i])
    Y[i] = np.array(Y[i])
#
plt.plot(X[0], Y[0], 'r')
plt.title(files[0])
plt.xscale('log')
plt.yscale('log')
plt.grid()
# plt.show()


# E = np.linspace(0.1, 10**4, 1000000) #MeV
OS = Cstop(X[0][45:], 8, 1, AO, IO)
CaS = Cstop(X[0][45:], 20, 1, ACa, ICa)
SS = Cstop(X[0][45:],16 , 1, AS, IS)
LiS = Cstop(X[0][45:], 3, 1, ALi, ILi)
BS = Cstop(X[0][45:], 5, 1, AB, IB)

plt.plot(X[0][45:], OS, 'go')
plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.show()

CaSO4 = fCa*CaS + fS*SS + fO4*OS
LTB = fLi*LiS + fB*BS + fO7*OS

plt.plot(X[0][45:], CaSO4)
plt.plot(X[0][45:], LTB)
plt.legend(['CaSO4', 'LTB'])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('MeV')
plt.ylabel('MeV cm^2/g')
plt.title('Stopping power')
plt.grid()
plt.show()

p = np.array([9, 2])

OS = Cstop(p, 8, 1, AO, IO)
CaS = Cstop(p, 20, 1, ACa, ICa)
SS = Cstop(p,16 , 1, AS, IS)
LiS = Cstop(p, 3, 1, ALi, ILi)
BS = Cstop(p, 5, 1, AB, IB)
CaSO4 = fCa*CaS + fS*SS + fO4*OS
LTB = fLi*LiS + fB*BS + fO7*OS
print(CaSO4)
print(LTB)
