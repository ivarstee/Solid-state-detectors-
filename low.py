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
    df = df.dropna(how='all').dropna(axis=1, how='all')
    array = df.to_numpy()
    x = np.zeros(616)
    y = np.zeros(616) #prøver dette så vi slipper å fikse alle data filer manuelt, 662 slutter x ved 350 grader celcius. Her burde vi slutte litt før
    for i in range(len(x)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    # n = np.where(y==np.nanmax(y))
    n = np.where(y[:400]==np.nanmax(y[:400])) 
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
    if N[i] > numb: 

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


files = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\low dose\*')

X = []
Y = []
N = []
for i in range(len(files)):
    x, y, n = reader(files[i])
    X.append(x)
    Y.append(y)
    N.append(n)
    # plt.plot(x, y)
    # plt.show()

Xmix = X[:36]
Ymix = Y[:36]
Nmix = N[:36]

XTLD100 = X[36:]
YTLD100 = Y[36:]

for i in range(len(Xmix)):
    index(Nmix, Ymix, i, 344)
    # plt.plot(Xmix[i], Ymix[i])
# plt.show()

ratio = np.zeros(len(Ymix))
Ca = np.zeros(len(Ymix))
LTB = np.zeros(len(Ymix))
for i in range(len(Ymix)):
    Ca[i] = np.nanmax(Ymix[i][:400])
    LTB[i] = np.nanmax(Ymix[i][400:556])
    n = np.where(Ymix[i][400:556]==np.nanmax(Ymix[i][400:556]))
    n = int(n[0]) + 400
    ratio[i] = np.nanmax(Ymix[i][:400])/np.nanmax(Ymix[i][400:556])
    print(files[i])
    plt.plot(Xmix[i], Ymix[i])
    plt.plot(Xmix[i][334], Ca[i], 'ro')
    plt.plot(Xmix[i][n], LTB[i], 'ro')
    # plt.title(files[i])
    plt.grid()
    plt.xlabel('Temperature [C]')
    plt.ylabel('intensity')
    plt.show()

D = np.array([100]*5 + [150]*6 + [25]*5 + [35]*5 + [45]*5 + [50]*5 + [75]*5)
D2 = np.array([100]*4 + [150]*4 + [25]*4 + [35]*4 + [45]*4 + [50]*4 + [75]*4)
TLD100max = []
for i in range(len(YTLD100)):
    TLD100max.append(np.nanmax(YTLD100[i]))

plt.plot(D, Ca, 'o')
plt.plot(D, LTB, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.ylabel('intensity')
plt.xlabel('Dose [mGy]')
plt.grid()
plt.show()

plt.plot(D, ratio, 'o') #vi har en usikkerhet i punktene på 25mGy siden vi egentlig ikke observerer noen LTB peak her
plt.ylabel('Ratio')
plt.xlabel('Dose [mGy]')
plt.grid()
plt.show()

plt.plot(D2, TLD100max, 'o')
plt.legend(['TLD100'])
plt.ylabel('intensity')
plt.xlabel('Dose [mGy]')
plt.grid()
plt.show()
