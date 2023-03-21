import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
import pandas as pd
import glob
from sklearn import linear_model
import scipy.stats



def reader(path):
    df = pd.read_excel(path, sheet_name='Sheet1')
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
    R = 1 - np.sum((y - fit)**2/np.sum((y - np.mean(y))**2))
    return R

def conf_plot(x, y, type):
    p_x = np.arange(np.min(x), np.max(x) + 1, 1)

    if type==1:
        alpha, beta = np.polyfit(x, y, type)
        print(alpha, beta)
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
        print(u)
        print(a, b)
        p = 2
        p_y1 = a*x2
        p_y = np.zeros(len(p_y1))

        for i in range(len(p_y1)):
            p_y[i] = np.sum(p_y1[i]) + b
            print(p_y[i], y[i])

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
    t = scipy.stats.t.ppf(q=1-0.025, df=len(x) - p)

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
    plow = p_y2 - abs(pred)
    phigh = p_y2 + abs(pred)



    plt.xlabel('Dose [Gy]')
    plt.ylabel('Amplitude')

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
    return [p_x, p_y2, lower, upper, plow, phigh]

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

CuLTB = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\TLD proton\CuLTB\*') #p1 = [:24], p2 = [24:]
TmCaSO4 = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\TLD proton\TmCaSO4\*') #p1 = [:20], p2 = [20:]

x_Cu_p1 = np.array([0.43]*3 + [0.61]*3 + [0.88]*3 + [1.18]*3 + [1.44]*3 + [2.10]*3 + [2.42]*3 + [3.19]*3)
x_Ca_p1 = np.array([0.34]*3 + [0.56]*3 + [0.90]*3 + [1.17]*2 + [1.63]*2 + [1.97]*3 + [2.39]*2 + [3.00]*2)
x_Cu_p2 = np.array([0.39]*3 + [0.64]*3 + [1.09]*3 + [1.31]*3 + [1.59]*3 + [2.11]*3 + [2.66]*3 + [3.00]*3)
x_Ca_p2 = np.array([0.44]*2 + [0.66]*2 + [0.90]*2 + [1.10]*2 + [1.59]*2 + [2.17]*2 + [2.55]*2 + [3.20]*2)

Cu = []
X_Cu = []
N_Cu = []
for i in range(len(CuLTB)):
    x, y, n = reader(CuLTB[i])
    X_Cu.append(x)
    Cu.append(np.amax(y))
    N_Cu.append(n)
    plt.plot(x, y)
    plt.grid()
    plt.title('CuLTB')
plt.show()
Cu_p1 = Cu[:24]
Cu_p2 = Cu[24:]

Ca = []
X_Ca = []
N_Ca = []
for i in range(len(TmCaSO4)):
    x, y, n = reader(TmCaSO4[i])
    X_Ca.append(x)
    Ca.append(np.amax(y))
    N_Ca.append(n)
    plt.plot(x, y)
    plt.grid()
    plt.title('TmCaSO4')
plt.show()

Ca_p1 = Ca[:20]
Ca_p2 = Ca[20:]

plt.plot(x_Ca_p1, Ca_p1, 'o')
plt.plot(x_Cu_p1, Cu_p1, 'o')
plt.grid()
plt.title('p1')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(x_Ca_p2, Ca_p2, 'o')
plt.plot(x_Cu_p2, Cu_p2, 'o')
plt.grid()
plt.title('p2')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()


# conf_Cu_p1 = conf_plot(x_Cu_p1, Cu_p1, 1)
# conf_Ca_p1 = conf_plot(x_Ca_p1, Ca_p1, 1)
# conf_Cu_p2 = conf_plot(x_Cu_p2, Cu_p2, 1)
# conf_Ca_p2 = conf_plot(x_Ca_p2, Ca_p2, 1)

conf_plot2(x_Ca_p1, Ca_p1)
conf_plot2(x_Cu_p1, Cu_p1)
conf_plot2(x_Ca_p2, Ca_p2)
conf_plot2(x_Cu_p2, Cu_p2)

# plt.xlabel('Dose [Gy]')
# plt.ylabel('Amplitude')
#
# plt.plot(x_Cu_p1, Cu_p1,'bo',label='CuLTB')
# plt.plot(conf_Cu_p1[0], conf_Cu_p1[1], 'k')
# plt.plot(x_Ca_p1, Ca_p1,'ro',label='TmCaSO4')
# plt.grid()
# plt.plot(conf_Cu_p1[0],conf_Cu_p1[2],'b--',label='confidence interval (95%)')
# plt.plot(conf_Cu_p1[0],conf_Cu_p1[3],'b--')
# plt.plot(conf_Cu_p1[0],conf_Cu_p1[-2],'g--',label='prediction interval (95%)')
# plt.plot(conf_Cu_p1[0],conf_Cu_p1[-1],'g--')
# # configure legend
# plt.legend(loc=0)
# leg = plt.gca().get_legend()
# ltext = leg.get_texts()
# plt.setp(ltext, fontsize=10)
#
# plt.show()

D = np.linspace(0.5, 3, 10)
ca1 = 344024*D - 16203
ca2 = 370109*D - 4565
ltb1 = 367108*D + 43893
ltb2 = 209476*D - 2787

plt.plot(D, ca1, 'o')
plt.plot(D, ltb1, 'o')
plt.grid()
plt.title('p1')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(D, ca2, 'o')
plt.plot(D, ltb2, 'o')
plt.grid()
plt.title('p2')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(D, ca1/ltb1, 'o')
plt.grid()
plt.title('p1')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()

plt.plot(D, ca2/ ltb2, 'o')
plt.grid()
plt.title('p2')
plt.legend(['CaSO4', 'LTB'])
plt.xlabel('Dose [Gy]')
plt.ylabel('Amplitude')
plt.show()
