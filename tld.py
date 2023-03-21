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


# har gjort de tre rene kurvene indentiske før lager gjennomsnitt

def reader(path):
    df = pd.read_excel(path, sheet_name='Sheet1')
    array = df.to_numpy()
    df = df.dropna(how='all').dropna(axis=1, how='all')
    array = df.to_numpy()
    x = np.zeros(1633) #1633 for 31.08
    y = np.zeros(1633) #662 for 28.08
    for i in range(len(x)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    n = np.where(y==np.nanmax(y))
    print(n)
    n = float(np.average(n[0]))


    return x, y, n


def r_sqr(fit, y):
    R = 1 - np.sum((y - fit)**2/np.sum((y - np.mean(y))**2))
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
        alpha, beta = np.polyfit(x, y, type)
        print(alpha, beta, r_sqr(alpha*x + beta, y))
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

        u = np.polyfit(x,y,2)
        a = u[:2]
        b = u[-1]
        # print(u)
        # print(a, b)
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
    plt.ylabel('a1/a2')

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


# files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Ravi_TLD\Excel\*')
files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\26.08 TLD\*')

#plott sammen ser stort sett ok ut, litt bekymret for ca4, om denne er forskøvet for langt til venstre. Prøver videre med gjennomsnitt for nå.

X = []
Y = []
N = []
for i in range(len(files)):
    x, y, n = reader(files[i])
    X.append(x)
    Y.append(y)
    N.append(n)
    # index(N, Y, i, 372.7)

for i in range(4, 8):
    index(N, Y, i, 476.25)
for i in range(0, 4):
    index(N, Y, i, 346.8)
for i in range(8, len(Y)):
    index(N, Y, i, 346.8)


#
# index(N, Y, 0, 339)
# index(N, Y, 2, 339)
# index(N, Y, 3, 472)
# index(N, Y, 5, 472)
#
#
# for i in range(6):
#     plt.plot(X[i], Y[i])
# plt.xlabel('degree celcius')
# plt.ylabel('intensity')
# plt.legend(['Ca1', 'Ca2', 'Ca4', 'Cu2', 'Cu3', 'Cu4'])
# plt.grid()
# plt.show()
# plt.show()
#


avxCa = (X[0] + X[1] + X[2] + X[3])/4
avyCa = (Y[0] + Y[1] + Y[2] + Y[3])/4
avxCu = (X[4] + X[5] + X[6] + X[7])/4
avyCu = (Y[4] + Y[5] + Y[6] + Y[7])/4
# avyCa = avyCa[:-12]
# avyCa = np.append([avyCa[0]]*12, avyCa)
# avyCu = avyCu[:-50]
# avyCu = np.append([avyCu[0]]*50, avyCu)#171 indeks = 79.6 grader i forskjell mellom de to toppene



AV = []
for i in range(len(avyCa)):
    AV.append([avyCa[i], avyCu[i]])

plt.plot(avxCa, AV)
plt.legend(['Average CaSO4', 'Average LTB'])
plt.grid()
plt.xlabel('degree celcius')
plt.ylabel('intensity')
plt.show()

A1 = []
A2 = []
Ca = []
Li = []
for i in range(8, len(X)):
    index(N, Y, i, 351)
    Ca.append(Y[i][351])
    Li.append(Y[i][522])
    plt.plot(X[i], Y[i])
    plt.plot(X[i][351], Y[i][351], 'ro')
    plt.plot(X[i][522], Y[i][522], 'ro')
    plt.grid()
    plt.title(files[i])
    plt.show()
#
#     #gammel kode!!
#     # a1, a2, r2 = fit(AV, X[i], Y[i], avyCa, avyCu, files[i])
#     # print(a1, a2)
#     # print(r2)
#     # A1.append(a1)
#     # A2.append(a2)
#
qx = np.array([100]*3 + [160]*3 + [225]*3)
eff = np.array([45]*4 + [67.0]*4 + [97.4]*4)

ratio = np.array(Ca)/np.array(Li)
# plt.plot(qx, ratio, 'o')
# plt.ylabel('intensity ratio (CaSO4/LTB)')
# plt.xlabel('kv')
# plt.grid()
# plt.show()
#

# plt.plot(np.delete(eff, -1), np.delete(ratio, 8), 'o')
plt.plot(eff, ratio, 'o')
plt.ylabel('intensity ratio (CaSO4/LTB)')
plt.xlabel('keV')
plt.grid()
plt.show()

# conf_plot(eff, ratio, 2)
# conf_plot(eff, ratio, 1)
conf_plot(np.delete(eff, -1), np.delete(ratio, 8), 1)
conf_plot(np.delete(eff, -1), np.delete(ratio, 8), 2)


plt.plot(eff, Li, 'o')
plt.plot(eff, Ca, 'o')
plt.ylabel('intensity')
plt.xlabel('keV')
plt.grid()
plt.legend(['LTB', 'CaSO4'])
plt.show()
#
#
# #begynner å se på stråling nr.2 fra 30.08.22 hvor vi har ny konsentrasjon av stoffer og ny oppvarmnings rate
#
#

files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Ravi_TLD(2)\Excel\*')
# files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\31.08 TLD\*')

X = []
Y = []
N = []
for i in range(len(files)):
    x, y, n = reader(files[i])
    y = y[:1633]
    x = x[:1633]
    X.append(x)
    Y.append(y)
    N.append(n)
N[13] = 710 #manuell korreksjon pga. feil

index(N, Y, 1, 735)
index(N, Y, 2, 1083)

for i in range(4):
    plt.plot(X[i], Y[i])
plt.xlabel('degree celcius')
plt.ylabel('intensity')
plt.legend(['Ca1', 'Ca2', 'Cu1', 'Cu4'])
plt.grid()
plt.show()
plt.show()

avxCa = (X[0] + X[1])/2
avyCa = (Y[0] + Y[1])/2
avxCu = (X[2] + X[3])/2
avyCu = (Y[2] + Y[3])/2

avyCa = avyCa[22:]
avyCa = np.append(avyCa, [avyCa[-1]]*22)
avyCu = avyCu[:-39]
avyCu = np.append([avyCu[0]]*39, avyCu)

AV = []
for i in range(len(avyCa)):
    AV.append([avyCa[i], avyCu[i]])

plt.plot(avxCa, AV)
plt.legend(['Average CaSO4', 'Average LTB'])
plt.grid()
plt.xlabel('degree celcius')
plt.ylabel('intensity')
plt.show()


A1 = []
A2 = []
Ca = []
Li = []
for i in range(4, 16):
    index(N, Y, i, 713)
    Ca.append(Y[i][713])
    Li.append(Y[i][1122])


# ACa = np.delete(np.array(ACa), 4)
# ALi = np.delete(np.array(ALi), 4)
qx = np.array([100]*4 + [160]*4 + [225]*4)
qx = np.delete(qx, 4)
eff = np.array([45]*3 + [67.0]*3 + [97.4]*4)



ratio = np.array(Ca)/np.array(Li)
ratio = np.delete(ratio, 4)
ratio = np.delete(ratio, 2)

plt.plot(eff, ratio, 'o')
plt.ylabel('intensity ratio (CaSO4/LTB)')
plt.xlabel('keV')
plt.title('amplitudes')
plt.grid()
plt.show()

conf_plot(eff, ratio, 2)
conf_plot(eff, ratio, 1)
#
#
#
#
# plt.plot(eff, Li, 'o')
# plt.plot(eff, Ca, 'o')
# plt.ylabel('intensity')
# plt.xlabel('keV')
# plt.grid()
# plt.legend(['LTB', 'CaSO4'])
# plt.show()



#
#
# X = []
# Y = []
# N = []
# for i in range(len(files)):
#     x, y, n = reader(files[i])
#     y = y[:1633]
#     x = x[:1633]
#     X.append(x)
#     Y.append(y)
#     N.append(n)
# N[13] = 710 #manuell korreksjon pga. feil
#
# index(N, Y, 1, 735)
# print(X[1][735])
# index(N, Y, 2, 1083)
# print(X[2][1083])
#
# for i in range(4):
#     plt.plot(X[i], Y[i])
# plt.xlabel('degree celcius')
# plt.ylabel('intensity')
# plt.legend(['Ca1', 'Ca2', 'Cu1', 'Cu4'])
# plt.grid()
# plt.show()
# plt.show()
#
# avxCa = (X[0] + X[1])/2
# avyCa = (Y[0] + Y[1])/2
# avxCu = (X[2] + X[3])/2
# avyCu = (Y[2] + Y[3])/2
#
# avyCa = avyCa[22:]
# avyCa = np.append(avyCa, [avyCa[-1]]*22)
# avyCu = avyCu[:-39]
# avyCu = np.append([avyCu[0]]*39, avyCu)
#
# AV = []
# for i in range(len(avyCa)):
#     AV.append([avyCa[i], avyCu[i]])
#
# plt.plot(avxCa, AV)
# plt.legend(['Average CaSO4', 'Average LTB'])
# plt.grid()
# plt.xlabel('degree celcius')
# plt.ylabel('intensity')
# plt.show()
#
#
# A1 = []
# A2 = []
# Ca = []
# Li = []
# for i in range(4, 16):
#     index(N, Y, i, 713)
#     Ca.append(Y[i][713])
#     Li.append(Y[i][1122])
#
#     # plt.plot(X[i], Y[i])
#     # plt.plot(X[i][713], Y[i][713], 'ro')
#     # plt.plot(X[i][1122], Y[i][1122], 'ro')
#     # plt.grid()
#     # plt.title(files[i])
#     # plt.show()
#     # a1, a2, r2 = fit(AV, X[i], Y[i], avyCa, avyCu, files[i])
#     # print(a1, a2)
#     # print(r2)
#     # A1.append(a1)
#     # A2.append(a2)
#
# # A1 = np.delete(A1, 4) #spør om vi kan fjerne dette punktet som en outlier
# # A2 = np.delete(A2, 4)
#
# #prøver annen måte å linærkombinere de rene spekterne
# ACa = []
# ALi = []
# for i in range(4, 16):
#     d = {'Sample': Y[i], 'Average_Ca': avyCa, 'Average_Li': avyCu}
#     df = pd.DataFrame(data=d)
#     y, X = dmatrices('Sample ~ Average_Ca + Average_Li', data=df, return_type='dataframe')
#     mod = sm.OLS(y, X)
#     res = mod.fit()
#     ACa.append(res.params[1])
#     ALi.append(res.params[2])
#     # print(res.params[1], res.params[2])
#     # print(res.summary())
#
# # ACa = np.delete(np.array(ACa), 4)
# # ALi = np.delete(np.array(ALi), 4)
# qx = np.array([100]*4 + [160]*4 + [225]*4)
# # qx = np.delete(qx, 4)
# eff = np.array([45]*4 + [67.0]*4 + [97.4]*4)
#
#
#
# ratio = np.array(Ca)/np.array(Li)
# # ratio = np.delete(ratio, 4)
#
# plt.plot(eff, ratio, 'o')
# plt.ylabel('intensity ratio (CaSO4/LTB)')
# plt.xlabel('keV')
# plt.title('amplitudes')
# plt.grid()
# plt.show()
#
#
# # plt.plot(qx, ACa/ALi, 'o')
# # plt.ylabel('intensity ratio (Ca/Li)')
# # plt.xlabel('kv')
# # plt.title('Fit')
# # plt.grid()
# # plt.show()
#
# # conf_plot(qx, ratio, 2)
# # conf_plot(qx, ratio, 1)
# x = sm.add_constant(eff)
# model =sm.OLS(ratio, x)
# res = model.fit()
# print(res.summary())
#
# #prøver å plotte regressjon med stats model + intervaller.
# prstd, iv_l, iv_u = wls_prediction_std(res)
# st, data, ss2 = summary_table(res, alpha=0.05)
#
# fittedvalues = data[:, 2]
# predict_mean_se  = data[:, 3]
# predict_mean_ci_low, predict_mean_ci_upp = data[:, 4:6].T
# predict_ci_low, predict_ci_upp = data[:, 6:8].T
#
# # Check we got the right things
# print(np.max(np.abs(res.fittedvalues - fittedvalues)))
# print(np.max(np.abs(iv_l - predict_ci_low)))
# print(np.max(np.abs(iv_u - predict_ci_upp)))
#
# # plt.plot(np.delete(eff, 4), ratio, 'bo', label='Data')
# # plt.plot(np.delete(eff, 4), fittedvalues, 'r-', lw=1.5, label='Regression line')
# # plt.plot(np.delete(eff, 4), predict_ci_low, 'g--', lw=1.5)
# # plt.plot(np.delete(eff, 4), predict_ci_upp, 'g--', lw=1.5, label='Prediction interval (95%)')
# # plt.plot(np.delete(eff, 4), predict_mean_ci_low, 'b--', lw=1.5)
# # plt.plot(np.delete(eff, 4), predict_mean_ci_upp, 'b--', label='Confidence interval (95%)')
# plt.plot(eff, ratio, 'bo', label='Data')
# plt.plot(eff, fittedvalues, 'r-', lw=1.5, label='Regression line')
# plt.plot(eff, predict_ci_low, 'g--', lw=1.5)
# plt.plot(eff, predict_ci_upp, 'g--', lw=1.5, label='Prediction interval (95%)')
# plt.plot(eff, predict_mean_ci_low, 'b--', lw=1.5)
# plt.plot(eff, predict_mean_ci_upp, 'b--', label='Confidence interval (95%)')
# plt.grid()
# plt.legend(loc=0)
# leg = plt.gca().get_legend()
# ltext = leg.get_texts()
# plt.setp(ltext, fontsize=10)
# plt.show()
#
# #intervallene ser identiske ut til de som jeg selv har programert
#
#
# plt.plot(eff, Li, 'o')
# plt.plot(eff, Ca, 'o')
# plt.ylabel('intensity')
# plt.xlabel('keV')
# plt.grid()
# plt.legend(['LTB', 'CaSO4'])
# plt.show()
#
# # plt.plot(qx, ACa, 'o')
# # plt.plot(qx, ALi, 'o')
# # plt.ylabel('intensity')
# # plt.xlabel('kv')
# # plt.grid()
# # plt.legend(['Li', 'Ca'])
# # plt.show()
# #Den nye metoden å linærkombinere de rene spekterne gir samme problem som tidligere at begge koeffisientene forandrer seg (ingen er konstante), trorlig pga. overlapp mellom kurvene.
# conf_plot(eff, ratio, 2)
# conf_plot(eff, ratio, 1)
