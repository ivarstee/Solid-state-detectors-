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
    y = np.zeros(662) #prøver dette så vi slipper å fikse alle data filer manuelt
    for i in range(len(x)):
        x[i] = array[i][0]
        y[i] = array[i][1]

    # n = np.where(y==np.nanmax(y))
    n = np.where(y[400:]==np.nanmax(y[400:])) #var nødvendig å trikse med indekser for at det skal funke i dette tilfellet.
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
            # print(p_y[i], y[i])

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



fade = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Gamma and proton\Fade\*')
proton = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Gamma and proton\Proton\*')
gamma = glob.glob(r'C:\Users\Hp\Desktop\skole\msc\Gamma and proton\Gamma\*')


P_X = []
P_Y = []
P_N = []
P_max = []
for i in range(len(proton)):
    x, y, n = reader(proton[i])
    plt.plot(x, y)
    plt.title(proton[i])
    plt.grid()
    plt.xlabel('Degrees [C]')
    plt.ylabel('intensity')
    plt.show()
    P_X.append(x)
    P_Y.append(y)
    P_N.append(n)
    P_max.append(np.amax(y))
P_max.pop(42) #fjerner en som ser ut som en outlier siden den ligger på utsiden av prediction interval

Ca_DP1 = np.array([0.59]*3 + [1.07]*3 + [1.46]*3 + [2.01]*3 + [3.00]*3)
Ca_DP2 = np.array([0.79]*3 + [1.06]*3 + [1.48]*3 + [1.92]*3 + [2.87]*3)
LTB_DP1 = np.array([0.47]*3 + [1.05]*3 + [1.48]*3 + [1.99]*3 + [2.98]*3)
LTB_DP2 = np.array([0.64]*3 + [1.09]*3 + [1.58]*3 + [2.15]*3 + [2.95]*3)

plt.plot(Ca_DP1, P_max[:15], 'o')
plt.plot(LTB_DP1, P_max[15:30], 'o')
plt.xlabel('Dose [Gy]')
plt.ylabel('intensity')
plt.title('P1')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.show()

plt.plot(Ca_DP2, P_max[30:45], 'o')
plt.plot(LTB_DP2, P_max[45:], 'o')
plt.xlabel('Dose [Gy]')
plt.ylabel('intensity')
plt.title('P2')
plt.legend(['CaSO4', 'LTB'])
plt.grid()
plt.show()


#conf_plot2 er mer ryddig men er bare linær. Confplot har mulighet til kvadratisk også ved å sette type = 2
conf_plot2(Ca_DP1, P_max[:15])
conf_plot2(LTB_DP1, P_max[15:30])
conf_plot2(Ca_DP2, P_max[30:45])
conf_plot2(LTB_DP2, P_max[45:])




#gamma data
D_gamma = np.array([0.5]*5 + [1.0]*5 + [1.5]*5 + [2.0]*5 + [2.5]*5 + [3.0]*5)
MI = [5, 6, 7, 8, 9, 15, 16, 17, 18, 19, 35, 36, 37, 38, 39, 50, 51, 52, 53, 54, 65, 66, 67, 68, 69, 80, 81, 82, 83, 84]
mixed = []

for i in MI:
    mixed.append(gamma[i])

del gamma[80:85]
del gamma[65:70]
del gamma[50:55]
del gamma[35:40]
del gamma[15:20]
del gamma[5:10]

#analyse av blandete prøver

mix_X = []
mix_Y = []
mix_N = []
for i in range(len(mixed)):
    x, y, n = reader(mixed[i])
    mix_X.append(x)
    mix_Y.append(y)
    mix_N.append(n)
    plt.plot(x, y)
    index(mix_N, mix_Y, i, 501) #indekser 12 og 19 må fikses manuelt siden toppene her har forskjellig størrelse forhold enn de andre kurvene.
plt.show()


ratio = np.zeros(len(mix_Y))
Ca = np.zeros(len(mix_Y))
LTB = np.zeros(len(mix_Y))
for i in range(len(mix_Y)):
    Ca[i] = np.nanmax(mix_Y[i][:400])
    LTB[i] = np.nanmax(mix_Y[i][400:])
    ratio[i] = np.nanmax(mix_Y[i][:400])/np.nanmax(mix_Y[i][400:])

plt.plot(D_gamma, ratio, 'o')
plt.xlabel('Dose [Gy]')
plt.ylabel('Ratio')
plt.legend(['CaSO4/LTB'])
plt.title('Mixed samples')
plt.grid()
plt.show()
conf_plot2(D_gamma, ratio)
plt.plot(D_gamma, Ca, 'o')
plt.plot(D_gamma, LTB, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.title('Mixed samples')
plt.xlabel('Dose [Gy]')
plt.ylabel('intensity')
plt.grid()
plt.show()

conf_plot2(D_gamma, Ca)
conf_plot2(D_gamma, LTB)


# for i in range(len(mix_Y)):
#     ind = np.where(y[:400]==np.nanmax(y[:400]))
#     ind = int(ind[0])
#     plt.plot(mix_X[i], mix_Y[i])
#     plt.plot(mix_X[i][ind], mix_Y[i][ind], 'ro')
#     # plt.plot(mix_X[i][501], mix_Y[i][501], 'ro')
#     plt.plot(mix_X[i][400], mix_Y[i][400], 'go')
#     plt.show()
#analyse av rene prøver

P_LTB = []

P_LTB.insert(0, gamma[-5:])
del gamma[-5:]
P_LTB.insert(0, gamma[-10:-5])
del gamma[-10:-5]
P_LTB.insert(0, gamma[-15:-10])
del gamma[-15:-10]
P_LTB.insert(0, gamma[-20:-15])
del gamma[-20:-15]
P_LTB.insert(0, gamma[-25:-20])
del gamma[-25:-20]
P_LTB.insert(0, gamma[-30:-25])
del gamma[-30:-25]

P_LTB = np.concatenate(P_LTB)
P_CaSO4 = gamma

max_LTB = []
max_CaSO4 = []

CaX = []
CaY = []
LTBX = []
LTBY = []

for i in range(len(P_LTB)):
    x1, y1, n1 = reader(P_CaSO4[i])
    x2, y2, n2 = reader(P_LTB[i])
    CaX.append(x1)
    CaY.append(y1)
    LTBX.append(x2)
    LTBY.append(y2)
    max_CaSO4.append(np.nanmax(y1))
    max_LTB.append(np.nanmax(y2))
    # plt.plot(x1, y1)
    # plt.plot(x2, y2)
    # plt.show()

max_LTB = np.array(max_LTB)
max_CaSO4 = np.array(max_CaSO4)

plt.plot(D_gamma, max_CaSO4/max_LTB, 'o')
plt.xlabel('Dose [Gy]')
plt.ylabel('Ratio')
plt.legend(['CaSO4/LTB'])
plt.title('Pure samples')
plt.grid()
plt.show()
conf_plot2(D_gamma, max_CaSO4/max_LTB)
plt.plot(D_gamma, max_CaSO4, 'o')
plt.plot(D_gamma, max_LTB, 'o')
plt.legend(['CaSO4', 'LTB'])
plt.title('Pure samples')
plt.xlabel('Dose [Gy]')
plt.ylabel('intensity')
plt.grid()
plt.show()

conf_plot2(D_gamma, max_CaSO4)
conf_plot2(D_gamma, max_LTB)

#start working on the fading results

"""
F_CaSO4 = fade[:10]
F_mix = fade[10:23]
F_LTB = fade[23:]
Fmax_CaSO4 = []
Fmax_LTB = []
Fmax_mix = []
for i in range(len(F_LTB)):
    x1, y1, n1 = reader(F_CaSO4[i])
    x2, y2, n2 = reader(F_LTB[i])
    Fmax_CaSO4.append(np.nanmax(y1))
    Fmax_LTB.append(np.nanmax(y2))
    plt.plot(x1, y1)
    plt.plot(CaX[i], CaY[i])
    plt.legend(['fade', 'original'])
    plt.title(F_CaSO4[i])
    plt.grid()
    plt.show()
    plt.plot(x2, y2)
    plt.plot(LTBX[i], LTBY[i])
    plt.legend(['fade', 'original'])
    plt.title(F_LTB[i])
    plt.grid()
    plt.show()

for i in range(len(F_mix)):
    x1, y1, n1 = reader(F_mix[i])
    # print(n1)
    # plt.plot(x1, y1)
    # plt.show()
    Fmax_mix.append(np.nanmax(y1))

fade_av = []
fade_av.append(np.average(Fmax_CaSO4[:5]))
fade_av.append(np.average(Fmax_CaSO4[5:]))
fade_av.append(np.average(Fmax_LTB[:5]))
fade_av.append(np.average(Fmax_LTB[5:]))
# fade_av.append(np.average(Fmax_mix[:5]))
# fade_av.append(np.average(Fmax_mix[5:10]))
# fade_av.append(np.average(Fmax_mix[10:]))

original_av = []
original_av.append(np.average(max_CaSO4[:5]))
original_av.append(np.average(max_CaSO4[5:10]))
original_av.append(np.average(max_LTB[:5]))
original_av.append(np.average(max_LTB[5:10]))

fade_av = np.array(fade_av)
original_av = np.array(original_av)

print(1 - fade_av/original_av)
"""
