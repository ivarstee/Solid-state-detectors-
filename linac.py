import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
import pandas as pd
import glob
from sklearn import linear_model
import scipy.stats
import statsmodels.api as sm
from patsy import dmatrices
from statsmodels.stats.outliers_influence import summary_table
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from scipy.interpolate import interp1d


files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\*.txt')
files = files[-2:]
files = files[::-1]
element = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\NIST\*")


X = [[], []]
Y = [[], []]
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
#bruker vektet gjennomsnitt for å finne de effektive energiene
eff = np.zeros(2)
for i in range(len(Y)):
    eff[i] = np.sum(X[i]*Y[i])/np.sum(Y[i])

muB = F[0](eff)
muCa = F[1](eff)
muH = F[2](eff)
muLi = F[3](eff)
muO = F[4](eff)
muS = F[5](eff)
muw = F[6](eff)

mu_CaSO4 = muCa*fCa + muS*fS + muO*fO4
mu_LTB = muLi*fLi + muB*fB + muO*fO7
comp_w = (2/10)*muH + (8/10)*muO

plt.plot(eff, mu_CaSO4/muw, 'o')
plt.plot(eff, mu_LTB/muw, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.xlabel('MeV')
plt.ylabel('(mu_en/ro)/(mu_en/ro)_w')
plt.show()
plt.plot(eff, mu_CaSO4/mu_LTB, 'o')
plt.grid()
plt.xlabel('MeV')
plt.ylabel('Ratio')
plt.show()

def spec_int(x, y, mu):
    D = 0
    dw = 0
    for i in range(len(y) - 1):
        D += y[i + 1]*(mu[i])*(x[i + 1] - x[i])
        dw += y[i + 1]*(F[6](x[i]))*(x[i + 1] - x[i])
    return D/dw


D_CaSO4 = np.zeros(2)
D_LTB = np.zeros(2)

for i in range(len(Y)):
    D_CaSO4[i] = spec_int(X[i]*10**(-3), Y[i], F[1](X[i]*10**(-3))*fCa + F[5](X[i]*10**(-3))*fS + F[4](X[i]*10**(-3))*fO4)
    D_LTB[i] = spec_int(X[i]*10**(-3), Y[i], F[3](X[i]*10**(-3))*fLi + F[0](X[i]*10**(-3))*fB + F[4](X[i]*10**(-3))*fO7)



plt.plot(eff, mu_CaSO4/mu_LTB, 'o')
plt.plot(eff, D_CaSO4/D_LTB, 'o')
plt.legend(['(mu_en/ro) ratio', 'Spectral average'])
plt.grid()
plt.xlabel('MeV')
plt.ylabel('Ratio')
plt.title('Ratio from x - ray spectra')
plt.show()

plt.plot(eff, D_CaSO4, 'o')
plt.plot(eff, D_LTB, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.title('Spektral integrasjon')
plt.xlabel('MeV')
# plt.title('Ratio from x - ray spectra')
plt.show()

def conf_plot2(x, y):
    p_x = np.arange(np.min(x), np.max(x) + 1, 1)
    p, V = np.polyfit(x, y, 1, cov="True")
    alpha, beta = p
    print("alpha: {} +/- {}".format(p[0], np.sqrt(V[0][0])))
    print("beta: {} +/- {}".format(p[1], np.sqrt(V[1][1])))
    p = 1
    p_y = alpha*x + beta
    p_y2 = alpha*p_x + beta
    plt.plot(p_x,alpha*p_x + beta,'r-',label='Regression line')
    plt.title('Linear regression and confidence limits')

    n = len(x)
    y_err = y - p_y
    RSS = np.sum(y_err**2)
    syx = np.sqrt(RSS/(n - p - 1))
    t = scipy.stats.t.ppf(q=1-0.025, df=len(x) - p - 1)
    print(t)
    mean_x = np.mean(x)
    # print(mean_x)
    confs = t*syx*np.sqrt(1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))
    pred = t*syx*np.sqrt(1 + 1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))


    lower = p_y2 - abs(confs)
    upper = p_y2 + abs(confs)

    plt.xlabel('kv')
    plt.ylabel('intensity')

    plt.plot(x,y,'bo',label='Sample observations')
    plt.grid()
    plt.plot(p_x,lower,'b--',label='confidence interval (95%)')
    plt.plot(p_x,upper,'b--')
    plt.plot(p_x,p_y2 - abs(pred),'g--',label='prediction interval (95%)')
    plt.plot(p_x,p_y2 + abs(pred),'g--')
    # configure legend
    plt.legend(loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=10)

    plt.show()

# conf_plot2(eff, mu_CaSO4/mu_LTB)
# conf_plot2(eff, D_CaSO4/D_LTB)
