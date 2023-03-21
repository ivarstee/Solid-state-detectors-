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


#ser på tld data fra 14.09.22

def reader(path):
    df = pd.read_excel(path, sheet_name='Sheet1')
    array = df.to_numpy()
    df = df.dropna(how='all').dropna(axis=1, how='all')
    array = df.to_numpy()
    x = np.zeros(len(array))
    y = np.zeros(len(array))
    for i in range(len(array)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    n = np.where(y==np.nanmax(y))
    n = float(n[0])


    return x, y, n


def r_sqr(fit, y):
    R = 1 - np.sum((y - fit)**2)/np.sum((y - np.mean(y))**2)
    return R

def fit(pure, x, y, sub1, sub2, file): #sub1 = avyCa, sub2 = avyCu
    reg = linear_model.LinearRegression().fit(pure, y)
    a1 = reg.coef_[0]
    a2 = reg.coef_[1]

    plt.plot(x, a1*sub1 + a2*sub2)
    plt.plot(x, y)
    plt.legend(['fit', 'data'])
    plt.grid()
    plt.xlabel('degree celcius')
    plt.ylabel('intensity')
    plt.title(file)
    plt.show()

    r2 = r_sqr(a1*sub1 + a2*sub2, y)
    return a1, a2, r2

def index(N, Y, i, numb):
    if N[i] > numb: #351 er gjennomsnitt av indeksen til originale topp punkter (i første måling), definerer derfor dette som indeks og bytter alle til dette.

        p = N[i]-numb
        Y[i] = Y[i][int(p):]
        Y[i] = np.append(Y[i], [Y[i][-1]]*int(p))

    if N[i] < numb:
        p = numb - N[i]
        Y[i] = Y[i][:-int(p)]
        Y[i] = np.append([Y[i][0]]*int(p), Y[i])

def conf_plot(x, y, type):
    p_x = np.arange(np.min(x), np.max(x) + 1, 1)

    if type==1:
        p_x = np.arange(np.min(x), np.max(x) + 1, 1)
        p, V = np.polyfit(x, y, 1, cov="True")
        alpha, beta = p
        print("alpha: {} +/- {}".format(p[0], np.sqrt(V[0][0])))
        print("beta: {} +/- {}".format(p[1], np.sqrt(V[1][1])))
        print(r_sqr(alpha*x + beta, y))
        print(1 - (((1 - r_sqr(alpha*x + beta, y))*(len(x) - 1))/(len(x) - 1 - 1)))
        p = 1
        p_y = alpha*x + beta
        p_y2 = alpha*p_x + beta
        plt.plot(p_x,alpha*p_x + beta,'r-',label='Regression line')
        plt.title('Linear regression and confidence limits')
    if type==2:
        p_x2 = []
        x2 = []
        for i in range(len(p_x)):
            p_x2.append(np.array([p_x[i]**2, p_x[i]]))
        for i in range(len(x)):
            x2.append(np.array([x[i]**2, x[i]]))

        u, V = np.polyfit(x,y,2, cov='True')
        a = u[:2]
        b = u[-1]
        # print(u)
        print(a, b)
        print("alpha: {} +/- {}".format(u[0], np.sqrt(V[0][0])))
        print("beta: {} +/- {}".format(u[1], np.sqrt(V[1][1])))
        print("gamma: {} +/- {}".format(u[2], np.sqrt(V[2][2])))
        p = 2
        p_y1 = a*x2
        p_y = np.zeros(len(p_y1))

        for i in range(len(p_y1)):
            p_y[i] = np.sum(p_y1[i]) + b
            # print(p_y[i], y[i])
        print(r_sqr(p_y, y))
        print(1 - (((1 - r_sqr(p_y, y))*(len(x) - 1))/(len(x) - 2 - 1)))
        p_y12 = a*p_x2
        p_y2 = np.zeros(len(p_y12))
        for i in range(len(p_y12)):
            p_y2[i] = np.sum(p_y12[i]) + b

        plt.plot(p_x,p_y2,'r-',label='Regression line')
        plt.title('quadratic regression and confidence limits')

    n = len(x)
    y_err = y - p_y
    RSS = np.sum(y_err**2)
    syx = np.sqrt(RSS/(n - p - 1))
    # mean_x = np.mean(x)
    t = scipy.stats.t.ppf(q=1-0.025, df=len(x) - p - 1)

    if type==1:
        mean_x = np.mean(x)
        # print(mean_x)
        confs = t*syx*np.sqrt(1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))
        pred = t*syx*np.sqrt(1 + 1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))
    if type==2:
        mean_x = np.mean(x2)
        confs1 = t*syx*np.sqrt(1.0/n + ((p_x2-mean_x)**2/np.sum((x2 - mean_x)**2)))
        pred1 = t*syx*np.sqrt(1 + 1.0/n + ((p_x2-mean_x)**2/np.sum((x2 - mean_x)**2)))
        confs = np.zeros(len(confs1))
        pred = np.zeros(len(confs1))
        for i in range(len(confs1)):
            confs[i] = np.sum(confs1[i])
            pred[i] = np.sum(pred1[i])
        # print(confs)


    lower = p_y2 - abs(confs)
    upper = p_y2 + abs(confs)


    plt.xlabel('keV')
    plt.ylabel('Ratio')

    plt.plot(x,y,'bo',label='Sample observations')
    plt.grid()
    plt.plot(p_x,lower,'b--',label='confidence interval (95%)')
    plt.plot(p_x,upper,'b--')
    plt.plot(p_x,p_y2 - abs(pred),'g--',label='prediction interval (95%)')
    plt.plot(p_x,p_y2 + abs(pred),'g--')
    # plt.errorbar(110.8, 1.6, yerr=0.3, xerr=15.6, fmt='ko', ecolor='k', elinewidth=2 ,capsize=10, capthick=2, label='Blind study result')
    # configure legend
    plt.legend(loc=0)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize=10)



    plt.show()
    # print(np.polyfit(p_x, p_y2 - abs(pred), 1),np.polyfit(p_x, p_y2 + abs(pred), 1))

files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\TLD (3)\*')

eff = np.array([45]*5 + [60.1]*5 + [75]*5 + [14.9]*5 + [105]*5 + [90.1]*5 + [120]*6 + [30]*5)

X = []
Y = []
N = []
for i in range(len(files)):
    x, y, n = reader(files[i])
    X.append(x)
    Y.append(y)
    N.append(n)
    index(N, Y, i, 334)


Ca = []
Li = []
ind = []

for i in range(len(Y)):
    Ca.append(Y[i][334])
    ind.append(int(np.where(Y[i][450:]==np.nanmax(Y[i][450:]))[0]) + 450)
    Li.append(Y[i][ind[i]])
    plt.plot(X[i], Y[i])
    plt.plot(X[i][334], Y[i][334], 'ro')
    plt.plot(X[i][ind[i]], Y[i][ind[i]], 'ro')
    plt.grid()
    plt.title(files[i])
    plt.show()

Li = np.array(Li)
Ca = np.array(Ca)
ratio = Ca/Li


plt.plot(eff, Li, 'o')
plt.plot(eff, Ca, 'o')
plt.xlabel('keV')
plt.ylabel('intensity')
plt.grid()
plt.legend(['LTB', 'CaSO4'])
plt.show()

plt.plot(np.delete(eff, 26), np.delete(ratio, 26), 'o')
plt.xlabel('keV')
plt.ylabel('intensity ratio (CaSO4/LTB)')
plt.grid()
plt.show()

HVL_Al = np.array([2.11]*5 + [5.83]*5 + [10.7]*5 + [0.243]*5 + [15.0]*5 + [12.7]*5 + [16.6]*6 + [1.22]*5) #mm
HVL_Cu = np.array([0.0698]*5 + [0.264]*5 + [0.730]*5 + [0.00814]*5 + [1.79]*5 + [1.17]*5 + [2.39]*6 + [0.0378]*5 ) #mm

plt.plot(HVL_Al, ratio, 'o')
plt.xlabel('HVL_AL [mm]')
plt.ylabel('intensity ratio (Ca/Li)')
plt.grid()
plt.show()
plt.plot(HVL_Cu, ratio, 'o')
plt.xlabel('HVL_Cu [mm]')
plt.ylabel('intensity ratio (Ca/Li)')
plt.grid()
plt.show()


conf_plot(np.delete(eff, 26), np.delete(ratio, 26), 1)
conf_plot(eff, ratio, 2)

x = sm.add_constant(eff)
model =sm.OLS(ratio, x)
res = model.fit()
# print(res.summary())

#prøver å plotte regressjon med stats model + intervaller.
prstd, iv_l, iv_u = wls_prediction_std(res)
st, data, ss2 = summary_table(res, alpha=0.05)

fittedvalues = data[:, 2]
predict_mean_se  = data[:, 3]
predict_mean_ci_low, predict_mean_ci_upp = data[:, 4:6].T
predict_ci_low, predict_ci_upp = data[:, 6:8].T

# Check we got the right things
# print(np.max(np.abs(res.fittedvalues - fittedvalues)))
# print(np.max(np.abs(iv_l - predict_ci_low)))
# print(np.max(np.abs(iv_u - predict_ci_upp)))
#
# plt.plot(eff, ratio, 'bo', label='Data')
# plt.plot(eff, fittedvalues, 'r-', lw=1, label='Regression line')
# plt.plot(eff, predict_ci_low, 'g--', lw=1)
# plt.plot(eff, predict_ci_upp, 'g--', lw=1, label='Prediction interval (95%)')
# plt.plot(eff, predict_mean_ci_low, 'b--', lw=1)
# plt.plot(eff, predict_mean_ci_upp, 'b--', lw=1,label='Confidence interval (95%)')
# plt.grid()
# plt.legend(loc=0)
# leg = plt.gca().get_legend()
# ltext = leg.get_texts()
# plt.setp(ltext, fontsize=10)
# plt.show()

# from scipy.interpolate import interp1d
# files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\*.txt')
# files = files[:-2]
# element = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\NIST\*")
# X = [[], [], [], [], [], [], [], []]
# Y = [[], [], [], [], [], [], [], []]
# int = np.zeros(len(files))
#
# for i in range(len(files)):
#     f = open(files[i], 'r')
#     for line in f:
#         row = line.split()
#         X[i].append(float(row[0]))
#         Y[i].append(float(row[1]))
#     X[i] = np.array(X[i])
#     Y[i] = np.array(Y[i])
#
#     plt.plot(X[i], Y[i])
#     plt.title(files[i])
#     plt.show()
#
# Ey = [[], [], [], [], [], [], []]
# Ex = [[], [], [], [], [], [], []]
# F = []
# for i in range(len(element)):
#     f = open(element[i], 'r')
#     for line in f:
#         row = line.split()
#         Ex[i].append(float(row[0]))
#         Ey[i].append(float(row[-1]))
#     Ex[i] = np.array(Ex[i])
#     Ey[i] = np.array(Ey[i])
#
#     F.append(interp1d(Ex[i], Ey[i], fill_value='extrapolate'))
#
#     plt.plot(Ex[i], F[i](Ex[i]), 'ro')
#     plt.plot(Ex[i], Ey[i])
#     plt.yscale('log')
#     plt.xscale('log')
#     plt.title(element[i])
#     plt.grid()
#     plt.show()
#
#
#
# fCa = 20/68
# fS = 16/68
# fO4 = 32/86
#
# fLi = 6/82
# fB = 20/82
# fO7 = 56/82
#
# # E = np.array([50, 60, 80, 100, 150])
# effT = np.array([14.9, 30, 45, 60.1, 75, 90.1, 105, 120])
# muB = F[0](effT*10**(-3))
# muCa = F[1](effT*10**(-3))
# muH = F[2](effT*10**(-3))
# muLi = F[3](effT*10**(-3))
# muO = F[4](effT*10**(-3))
# muS = F[5](effT*10**(-3))
# muw = F[6](effT*10**(-3))
#
#
#
# mu_CaSO4 = muCa*fCa + muS*fS + muO*fO4
# mu_LTB = muLi*fLi + muB*fB + muO*fO7
# comp_w = (2/10)*muH + (8/10)*muO

# plt.plot(eff, mu_CaSO4/muw, 'o')
# plt.plot(eff, mu_LTB/muw, 'o')
# plt.legend(['CaSO4', 'LTB'])
# plt.grid()
# plt.xlabel('keV')
# plt.ylabel('(mu_en/ro)/(mu_en/ro)_w')
# plt.title('mass energy absorption coefficients')
# plt.show()


# def spec_int(x, y, mu):
#     D = 0
#     dw = 0
#     for i in range(len(y) - 1):
#         D += y[i + 1]*(mu[i])*(x[i + 1] - x[i])
#         dw += y[i + 1]*(F[6](x[i]))*(x[i + 1] - x[i])
#     return D/dw
#
#
# D_CaSO4 = np.zeros(8)
# D_LTB = np.zeros(8)
#
# # print(spec_int(X[0], Y[0], F[1](X[0])*fCa + F[4](X[0])*fS + F[3](X[0])*fO4))
#
# for i in range(len(Y)):
#     D_CaSO4[i] = spec_int(X[i]*10**(-3), Y[i], F[1](X[i]*10**(-3))*fCa + F[5](X[i]*10**(-3))*fS + F[4](X[i]*10**(-3))*fO4)
#     D_LTB[i] = spec_int(X[i]*10**(-3), Y[i], F[3](X[i]*10**(-3))*fLi + F[0](X[i]*10**(-3))*fB + F[4](X[i]*10**(-3))*fO7)
#
#
# plt.plot(effT, mu_CaSO4/mu_LTB, 'o')
# plt.plot(effT, D_CaSO4/D_LTB, 'o')
# plt.plot(eff, ratio, 'o')
# plt.legend(['(mu_en/ro) ratio', 'Spectral average', 'Data'])
# plt.grid()
# plt.xlabel('keV')
# plt.ylabel('Ratio')
# plt.title('Ratio from x - ray spectra')
# plt.show()

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


new = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\10 kV\*')
X_new = []
Y_new = []
N_new = []
for i in range(len(new)):
    x, y, n = reader(new[i])
    x = x[:660]
    y = y[:660]
    X_new.append(x)
    Y_new.append(y)
    N_new.append(n)
    index(N_new, Y_new, i, 343)
    plt.plot(X_new[i], Y_new[i])
plt.show()

Ca_new = []
LTB_new = []

for i in range(len(X_new)):
    Ca_new.append(Y_new[i][343])
    LTB_new.append(np.nanmax(Y_new[i][470:570]))


plt.plot(eff, ratio, 'go')
plt.plot([7.01, 7.01, 7.01, 7.01, 7.01], (np.array(Ca_new)/np.array(LTB_new)), 'ro')
plt.xlabel('keV')
plt.ylabel('intensity ratio (CaSO4/LTB)')
plt.grid()
plt.show()

sam = np.array([np.average(np.array(Ca_new)/np.array(LTB_new))] + [np.average(ratio[15:20])] + [np.average(ratio[36:])] + [np.average(ratio[:5])] + [np.average(ratio[5:10])] + [np.average(ratio[10:15])] + [np.average(ratio[25:30])] + [np.average(ratio[20:25])]  + [np.average(ratio[30:36])])/2.54413344386335

E = np.array([7.01, 14.9, 30, 45, 60.1, 75, 90.1, 105, 120])
mono = -7.90507810e-08*E**4 + 2.62515418e-05*E**3 -2.91822638e-03*E**2 + 1.06353740e-01*E + 5.77044726e-01
av = -5.57834805e-08*E**4 + 1.89994981e-05*E**3 -2.27334706e-03*E**2 + 9.65888059e-02*E + 2.54765906e-01
plt.plot(E, sam, 'o')

# plt.plot(E, [0.67062341, 0.88243076, 1.        , 0.87019004, 0.59100114, 0.40386731, 0.26474418, 0.1916607 , 0.16211845])
# plt.plot(E, [0.52193247, 0.87308737, 1.        , 0.97907359, 0.87082139, 0.58748485, 0.50323908, 0.34853298, 0.24861723])
plt.plot(E, mono, 'r',label='Monoenergetic 4. deg. polynomial curve')
plt.plot(E, av, 'b',label='Spectral average 4. deg. polynomial curve')
plt.xlabel('keV')
plt.ylabel('Normalized intensity ratio (CaSO4/LTB)')
# plt.legend(['Monoenergetic 4. deg. polynomial curve', 'Spectral average 4. deg. polynomial curve'])
# plt.plot(eff, ratio/2.54413344386335, 'go')
# plt.plot([7.01, 7.01, 7.01, 7.01, 7.01], (np.array(Ca_new)/np.array(LTB_new))/2.54413344386335, 'go', label='Sample data')
plt.grid()
plt.legend(loc='lower left')
leg = plt.gca().get_legend()
ltext = leg.get_texts()
plt.setp(ltext, fontsize=10)
plt.show()

# print(r_sqr(sam, np.array([0.67062341, 0.88243076, 1.        , 0.87019004, 0.59100114, 0.40386731, 0.26474418, 0.1916607 , 0.16211845])))
# print(r_sqr(sam, np.array([0.52193247, 0.87308737, 1.        , 0.97907359, 0.87082139, 0.58748485, 0.50323908, 0.34853298, 0.24861723])))

print(r_sqr(mono, sam))
print(r_sqr(av, sam))

def radj(F, Y, d):
    R = r_sqr(F, Y)
    print(1 - (((1 - R)*(9 - 1))/(9 - d - 1)))

def quad_diff(y, fit):
    return np.sum((y - fit)**2)

print(quad_diff(mono, sam))
print(quad_diff(av, sam))
