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
    x = np.zeros(662)
    y = np.zeros(662)
    for i in range(len(x)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    # n = np.where(y==np.nanmax(y))
    n = np.where(y[400:]==np.nanmax(y[400:]))
    n = float(n[0] + 400)


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
    # print(t)
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


data = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Aarhus data\*')


max = []
for i in range(len(data)):
    x, y, n = reader(data[i])
    plt.plot(x, y)
    plt.title(data[i])
    plt.grid()
    plt.xlabel('Degrees [C]')
    plt.ylabel('intensity')
    plt.show()
    max.append(np.nanmax(y))

CaSO4 = max[:24]
LTB = max[24:]

for i in range(len(CaSO4)):
    x, y, n = reader(data[i])
    print(data[i])
    plt.plot(x, y)
    plt.grid()
    plt.xlabel('Degrees [C]')
    plt.ylabel('intensity')
plt.show()

for i in range(len(LTB)):
    x, y, n = reader(data[i + 24])
    print(data[i + 24])
    plt.plot(x, y)
    plt.grid()
    plt.xlabel('Degrees [C]')
    plt.ylabel('intensity')
plt.show()


Ca_pos = np.array([2]*6 + [3]*6 + [4]*6 + [1]*6)
LTB_pos = np.array([2]*6 + [3]*6 + [4]*6 + [1]*6)
lab = ['1a + 3a', '1b', '1c', '2a']
lab = np.array(['1a + 3a']*6 + ['1b']*6 + ['1c']*6 + ['2a']*6)

plt.plot(Ca_pos, CaSO4, 'o')
plt.plot(LTB_pos, LTB, 'o')
plt.xticks(Ca_pos, lab)
plt.grid()
plt.xlabel('Position')
plt.ylabel('intensity')
plt.legend(['CaSO4', 'LTB'])
plt.show()

plt.plot(LTB_pos, LTB, 'o')
plt.grid()
plt.xlabel('Position')
plt.ylabel('intensity')
plt.legend(['LTB'])
plt.show()

plt.plot(Ca_pos, np.array(CaSO4)/np.array(LTB), 'o')
plt.xticks(Ca_pos, lab)
plt.grid()
plt.xlabel('Position')
plt.ylabel('intensity ratio')
plt.legend(['CaSO4/LTB'])
plt.show()

av = []
av.append(np.average(np.array(CaSO4[:6])/np.array(LTB[:6])))
av.append(np.average(np.array(CaSO4[6:12])/np.array(LTB[6:12])))
av.append(np.average(np.array(CaSO4[12:18])/np.array(LTB[12:18])))
av.append(np.average(np.array(CaSO4[18:24])/np.array(LTB[18:24])))

plt.plot(np.array([2] + [3] + [4] + [1]), av, 'o')
plt.xticks(np.array([2] + [3] + [4] + [1]), np.array(['1a + 3a'] + ['1b'] + ['1c'] + ['2a']))
plt.grid()
# plt.ylim(0.8, 1.15)
plt.xlabel('Position')
plt.ylabel('intensity ratio')
plt.legend(['Average CaSO4/LTB'])
plt.show()
