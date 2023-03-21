import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
import pandas as pd
import glob
import scipy.stats
from sklearn import linear_model

t = [5, 5, 5, 10, 10, 10, 15, 15, 15, 20, 20, 20]
time = np.linspace(5, 20, 100)

# max = np.array([50, 70, 100, 120, 135, 150, 180, 200])
# eff = np.array([99.9, 86.7, 78.5, 67.2, 60.2, 51.6, 41.3, 34.3])
eff = np.array([34.3, 41.3, 51.6, 60.2, 67.2, 78.5, 86.7, 99.9]) #forsøk 1
Pu = [1.04, 1.04, 1.03, 1.03, 1.03, 1.02, 1.02, 1.02]
mu = [1.012, 1.017, 1.027, 1.036, 1.045, 1.058, 1.068, 1.075]

mu_coef = np.polyfit(eff, mu, 1)
Pu_coef = np.polyfit(eff, Pu, 1)
plt.plot(eff, Pu, 'o')
plt.plot(eff, Pu_coef[0]*eff + Pu_coef[1])
plt.show()
plt.plot(eff, mu, 'o')
plt.plot(eff, mu_coef[0]*eff + mu_coef[1])
plt.show()

N = 43.77
Ku = 1

# kV120 = np.array([0.21, 0.28, 0.28, 0.94, 0.92, 0.92, 1.58, 1.60, 1.54, 2.26, 2.31, 2.30, 7.60, 7.60, 7.61])
# kV105 = np.array([0.37, 0.29, 0.33, 1.20, 1.21, 1.16, 2.11, 2.04, 2.14, 3.05, 3.09, 2.97, 10.13, 10.27, 10.12])
# kV90 = np.array([0.54, 0.45, 0.64, 1.86, 1.97, 1.84, 3.35, 3.22, 3.33, 4.65, 4.77, 4.58, 15.65, 15.57, 15.69])
# kV75 = np.array([0.20, 0.15, 0.21, 0.62, 0.62, 0.61, 1.10, 1.03, 1.03, 1.54, 1.54, 1.52, 5.04, 5.00, 5.05])
# kV60 = np.array([0.44, 0.38, 0.47, 1.35, 1.32, 1.27, 2.22, 2.12, 2.09, 3.02, 2.89, 2.96, 9.82, 9.82, 9.86])
# kV45 = np.array([0.58, 0.59, 0.49, 1.56, 1.54, 1.52, 2.42, 2.46, 2.48, 3.49, 3.39, 3.33, 11.18, 11.02, 11.11])
# kV30 = np.array([0.16, 0.15, 0.18, 0.44, 0.42, 0.43, 0.74, 0.70, 0.74, 1.02, 1.03, 0.99, 3.27, 3.29, 3.30])
# kV15 = np.array([0.01, 0.02, 0.02, 0.07, 0.06, 0.07, 0.11, 0.11, 0.11, 0.16, 0.15, 0.16, 0.53, 0.53, 0.53])
#
#
#
def Dw(M, N, k, mu_ro, p):
    # Ktp = 1.031968 # til første forsøk med 23.1 grader
    # Ktp = 1.0366 ved std og CU dosimetri
    # Ktp = 1.0256 #ved 150 kV måling
    Ktp = 0.9809 #alanin dosimetri
    return M*N*k*mu_ro*p*Ktp

kV_STD = np.array([0.6, 0.56, 0.41, 1.84, 1.81, 1.68, 3.02, 3.04, 3.11, 4.17, 4.30, 4.28])
DSTD = Dw(kV_STD, N, Ku, mu_coef[0]*97.4 + mu_coef[1], Pu_coef[0]*97.4 + Pu_coef[1])
DSTD_coef = np.polyfit(t, DSTD, 1)
plt.plot(t, DSTD, 'o')
plt.plot(time, DSTD_coef[0]*time + DSTD_coef[1])
plt.title('97.4 keV')
plt.grid()
plt.show()

print(DSTD_coef[0],DSTD_coef[0]*60)
print(10/(DSTD_coef[0]))

#
# D120 = Dw(kV120, N, Ku, mu_coef[0]*120 + mu_coef[1], Pu_coef[0]*120 + Pu_coef[1])
# D120_coef = np.polyfit(t, D120[:-3], 1)
# plt.plot(t, D120[:-3], 'o')
# plt.plot(time, D120_coef[0]*time + D120_coef[1])
# plt.title('120 keV')
# plt.grid()
# plt.show()
#
# print(D120_coef[0],D120_coef[0]*60 ,np.average(D120[-3:]))
# print(10/(D120_coef[0]))
#
# D105 = Dw(kV105, N, Ku, mu_coef[0]*105 + mu_coef[1], Pu_coef[0]*105 + Pu_coef[1])
# D105_coef = np.polyfit(t, D105[:-3], 1)
# plt.plot(t, D105[:-3], 'o')
# plt.plot(time, D105_coef[0]*time + D105_coef[1])
# plt.title('105 keV')
# plt.grid()
# plt.show()
#
# print(D105_coef[0],D105_coef[0]*60 ,np.average(D105[-3:]))
# print(10/(D105_coef[0]))
#
# D90 = Dw(kV90, N, Ku, mu_coef[0]*90 + mu_coef[1], Pu_coef[0]*90 + Pu_coef[1])
# D90_coef = np.polyfit(t, D90[:-3], 1)
# plt.plot(t, D90[:-3], 'o')
# plt.plot(time, D90_coef[0]*time + D90_coef[1])
# plt.title('90 keV')
# plt.grid()
# plt.show()
#
# print(D90_coef[0],D90_coef[0]*60 ,np.average(D90[-3:]))
# print(10/(D90_coef[0]))
#
# D75 = Dw(kV75, N, Ku, mu_coef[0]*75 + mu_coef[1], Pu_coef[0]*75 + Pu_coef[1])
# D75_coef = np.polyfit(t, D75[:-3], 1)
# plt.plot(t, D75[:-3], 'o')
# plt.plot(time, D75_coef[0]*time + D75_coef[1])
# plt.title('75 keV')
# plt.grid()
# plt.show()
#
# print(D75_coef[0],D75_coef[0]*60 ,np.average(D75[-3:]))
# print(10/(D75_coef[0]))
#
# D60 = Dw(kV60, N, Ku, mu_coef[0]*60 + mu_coef[1], Pu_coef[0]*60 + Pu_coef[1])
# D60_coef = np.polyfit(t, D60[:-3], 1)
# plt.plot(t, D60[:-3], 'o')
# plt.plot(time, D60_coef[0]*time + D60_coef[1])
# plt.title('60 keV')
# plt.grid()
# plt.show()
#
# print(D60_coef[0],D60_coef[0]*60 ,np.average(D60[-3:]))
# print(10/(D60_coef[0]))
#
# D45 = Dw(kV45, N, Ku, mu_coef[0]*45 + mu_coef[1], Pu_coef[0]*45 + Pu_coef[1])
# D45_coef = np.polyfit(t, D45[:-3], 1)
# plt.plot(t, D45[:-3], 'o')
# plt.plot(time, D45_coef[0]*time + D45_coef[1])
# plt.title('45 keV')
# plt.grid()
# plt.show()
#
# print(D45_coef[0],D45_coef[0]*60 ,np.average(D45[-3:]))
# print(10/(D45_coef[0]))
#
# D30 = Dw(kV30, N, Ku, mu_coef[0]*30 + mu_coef[1], Pu_coef[0]*30 + Pu_coef[1])
# D30_coef = np.polyfit(t, D30[:-3], 1)
# plt.plot(t, D30[:-3], 'o')
# plt.plot(time, D30_coef[0]*time + D30_coef[1])
# plt.title('30 keV')
# plt.grid()
# plt.show()
#
# print(D30_coef[0],D30_coef[0]*60 ,np.average(D30[-3:]))
# print(10/(D30_coef[0]))
#
# D15 = Dw(kV15, N, Ku, mu_coef[0]*15 + mu_coef[1], Pu_coef[0]*15 + Pu_coef[1])
# D15_coef = np.polyfit(t, D15[:-3], 1)
# plt.plot(t, D15[:-3], 'o')
# plt.plot(time, D15_coef[0]*time + D15_coef[1])
# plt.title('15 keV')
# plt.grid()
# plt.show()
#
# print(D15_coef[0],D15_coef[0]*60 ,np.average(D15[-3:]))
# print(10/(D15_coef[0]))


#ny dosimetri til forsøk 2

max = np.array([100, 225])
doses = np.array([5, 10, 15, 25, 50, 75, 100])

Cu = np.array([0.01, 0.01, 0.02, 0.06, 0.06, 0.05, 0.10, 0.10, 0.10, 0.13, 0.14, 0.14, 0.46, 0.46, 0.47])
std = np.array([0.60, 0.56, 0.41, 1.84, 1.81, 1.68, 3.02, 3.04, 3.11, 4.17, 4.30, 4.28, 14.16, 14.35, 14.33])

DCu = Dw(Cu, N, Ku, mu_coef[0]*74.9 + mu_coef[1], Pu_coef[0]*74.9 + Pu_coef[1])
DCu_coef = np.polyfit(t, DCu[:-3], 1)
plt.plot(t, DCu[:-3], 'o')
plt.plot(time, DCu_coef[0]*time + DCu_coef[1])
plt.title('Cu set-up')
plt.grid()
plt.show()

print(DCu_coef[0],DCu_coef[0]*60 ,np.average(DCu[-3:]))
print(doses/(DCu_coef[0])*4)# ganger på 4 siden vi målte med 20 mA men bruker 5 mA i forsøket


Dstd = Dw(std, N, Ku, mu_coef[0]*97.4 + mu_coef[1], Pu_coef[0]*97.4 + Pu_coef[1])
Dstd_coef = np.polyfit(t, Dstd[:-3], 1)
plt.plot(t, Dstd[:-3], 'o')
plt.plot(time, Dstd_coef[0]*time + Dstd_coef[1])
plt.title('225 kV standard filtration')
plt.grid()
plt.show()

print(Dstd_coef[0],Dstd_coef[0]*60 ,np.average(Dstd[-3:]))
print((doses/(Dstd_coef[0]))*10) #ganger med 10 siden vi målte på 10 mA, men bruker 1 mA under forsøket

kV160 = np.array([0.59, 0.68, 0.72, 1.86, 1.70, 1.76, 3.03, 3.06, 2.94, 4.22, 4.06, 4.11, 13.83, 13.75, 13.65])

D160 = Dw(std, N, Ku, mu_coef[0]*67 + mu_coef[1], Pu_coef[0]*67 + Pu_coef[1])
D160_coef = np.polyfit(t, D160[:-3], 1)
plt.plot(t, D160[:-3], 'o')
plt.plot(time, D160_coef[0]*time + D160_coef[1])
plt.title('160 kV 2.02+0.1')
plt.grid()
plt.show()

print(D160_coef[0],D160_coef[0]*60 ,np.average(D160[-3:]))


#alanin dosimetri
D25 = np.array([6.82, 6.80, 6.73])
D48 = np.array([14.58, 14.42, 14.64])
D72 = np.array([13.59, 13.66, 13.56])
D96 = np.array([9.11, 9.13, 9.12])
D120 = np.array([7.92, 7.95, 7.96])

print(np.average(Dw(D25, N, Ku, mu_coef[0]*25.09 + mu_coef[1], Pu_coef[0]*25.09 + Pu_coef[1]))/1000)
print(np.average(Dw(D48, N, Ku, mu_coef[0]*48.8 + mu_coef[1], Pu_coef[0]*48.8 + Pu_coef[1]))/1000)
print(np.average(Dw(D72, N, Ku, mu_coef[0]*72.7 + mu_coef[1], Pu_coef[0]*72.7 + Pu_coef[1]))/1000)
print(np.average(Dw(D96, N, Ku, mu_coef[0]*96.6 + mu_coef[1], Pu_coef[0]*96.6 + Pu_coef[1]))/1000)
print(np.average(Dw(D120, N, Ku, mu_coef[0]*120 + mu_coef[1], Pu_coef[0]*120 + Pu_coef[1]))/1000)
