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


def reader(path):
    df = pd.read_excel(path, sheet_name='Sheet1')
    df = df.dropna(how='all').dropna(axis=1, how='all')
    array = df.to_numpy()
    x = np.zeros(662)
    y = np.zeros(662) #prøver dette så vi slipper å fikse alle data filer manuelt
    for i in range(len(x)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    # n = np.where(y==np.nanmax(y))
    n = np.where(y[:425]==np.nanmax(y[:425])) #var nødvendig å trikse med indekser for at det skal funke i dette tilfellet.
    n = float(n[0])


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

def conf_plot2(x, y):
    p_x = np.arange(np.min(x), np.max(x) + 1, 1)
    p, V = np.polyfit(x, y, 1, cov="True")
    alpha, beta = p
    print("alpha: {} +/- {}".format(p[0], np.sqrt(V[0][0])))
    print("beta: {} +/- {}".format(p[1], np.sqrt(V[1][1])))
    print(r_sqr(alpha*x + beta, y))
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
    mean_x = np.mean(x)
    # print(mean_x)
    confs = t*syx*np.sqrt(1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))
    pred = t*syx*np.sqrt(1 + 1.0/n + ((p_x-mean_x)**2/np.sum((x - mean_x)**2)))


    lower = p_y2 - abs(confs)
    upper = p_y2 + abs(confs)

    plt.xlabel('Dose')
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
    return p_x, lower, upper, p_y2 - abs(pred), p_y2 + abs(pred), alpha, beta

files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\linac\*')
TLD100 = []
TLD100.append(files[25:45])
TLD100.append(files[70:])
del files[70:]
del files[25:45]
TLD100 = np.concatenate(TLD100)

X = []
Y = []
N = []
for i in range(len(files)):
    x, y, n = reader(files[i])
    X.append(x)
    Y.append(y)
    N.append(n)
    index(N, Y, i, 356)



# her skalerer vi responsen for å bruke den til Blindstudie
ratio = np.zeros(len(Y))
Ca = np.zeros(len(Y))
LTB = np.zeros(len(Y))
for i in range(len(Y)):
    Ca[i] = np.nanmax(Y[i][:425])*3.384998722046466
    LTB[i] = np.nanmax(Y[i][425:])#*1.2379433800666064
    # ratio[i] = np.nanmax(Y[i][:425])/np.nanmax(Y[i][425:])


D = np.array([0.2]*5 + [0.4]*5 + [0.6]*5 + [0.8]*5 + [1.0]*5)

# plt.plot(D, ratio[:25], 'o')
# plt.xlabel('Dose [Gy]')
# plt.ylabel('Ratio')
# plt.legend(['CaSO4/LTB'])
# plt.title('15 MV')
# plt.grid()
# plt.show()

plt.plot(D, Ca[:25], 'o')
plt.plot(D, LTB[:25], 'o')
plt.xlabel('Dose [Gy]')
plt.ylabel('intensity')
plt.title('15 MV')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.show()

# print('15 MV ratio')
# conf_plot2(D, ratio[:25])
# print('15 MV CaSO4')
# conf_plot2(D, Ca[:25])
# print('15 MV LTB')
# conf_plot2(D, LTB[:25])

dummy = 400000
#
# def find_dose(Ca_peak, LTB_peak):
#     px, Ca_CL, Ca_CU, Ca_PL, Ca_PU, Ca_Malpha, Ca_Mbeta = conf_plot2(D, Ca[:25])
#     Ca_CLalpha, Ca_CLbeta = np.polyfit(px, Ca_CL, 1)
#     Ca_CUalpha, Ca_CUbeta = np.polyfit(px, Ca_CU, 1)
#     Ca_PLalpha, Ca_PLbeta = np.polyfit(px, Ca_PL, 1)
#     Ca_PUalpha, Ca_PUbeta = np.polyfit(px, Ca_PU, 1)
#
#
#     px, LTB_CL, LTB_CU, LTB_PL, LTB_PU, LTB_Malpha, LTB_Mbeta = conf_plot2(D, LTB[:25])
#     LTB_CLalpha, LTB_CLbeta = np.polyfit(px, LTB_CL, 1)
#     LTB_CUalpha, LTB_CUbeta = np.polyfit(px, LTB_CU, 1)
#     LTB_PLalpha, LTB_PLbeta = np.polyfit(px, LTB_PL, 1)
#     LTB_PUalpha, LTB_PUbeta = np.polyfit(px, LTB_PU, 1)
#     return [(Ca_peak-Ca_PLbeta)/Ca_PLalpha, (Ca_peak-Ca_PUbeta)/Ca_PUalpha, (Ca_peak-Ca_Mbeta)/Ca_Malpha], [(LTB_peak-LTB_PLbeta)/LTB_PLalpha, (LTB_peak-LTB_PUbeta)/LTB_PUalpha, (LTB_peak-LTB_Mbeta)/LTB_Malpha]


def find_energy(ratio):
    return (ratio - 3.49903486)/(-0.02341907), (ratio - 4.81288693)/(-0.02344136), (ratio - 4.155960895267299)/(-0.023430216453664908)

X = []
Y = []
N = []
Ca_sample = []
LTB_sample = []
files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Bilnd 1\*')
for i in range(len(files)):
    x, y, n = reader(files[i])
    X.append(x)
    Y.append(y)
    N.append(n)
    index(N, Y, i, 358)
    Ca_sample.append(Y[i][358])
    LTB_sample.append(np.nanmax(Y[i][460:]))

# print(find_dose(np.average(Ca_sample), np.average(LTB_sample)))
print(find_energy(np.average(Ca_sample)/np.average(LTB_sample)))
