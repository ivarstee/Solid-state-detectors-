import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import argrelextrema
import pandas as pd
import glob
import scipy.stats
from sklearn import linear_model


class finder:

    def gaussian(x, a, mu, Sigma):
        return a*np.exp(-np.power(x - mu, 2.)/(2*np.power(Sigma, 2.)))

    def der(x, y):
        d = np.zeros(len(y))
        for i in range(len(y) - 1):
            d[i] = (y[i + 1] - y[i])/(x[i + 1] - x[i])
        return d


    def bimodal_der(X, a1, mu1, a2, bias):
        x, u = X
        g = finder.gaussian(x, a1, mu1, 3.7377) + finder.gaussian(x, a2, mu1-8, 1.9795) # har prøvd -5 siden dette er avstanden mellom ekstremal punktene i 500 Gy 10 dB
        d = np.zeros(len(g))
        for i in range(len(g) - 1):
            d[i] = (g[i + 1] - g[i])/(x[i + 1] - x[i])

        return d + u + bias

    def r_sqr(fit, y):
        R = 1 - np.sum((y - fit)**2/np.sum((y - np.mean(y))**2))
        return R

    def smooth(y):
        y_r = y
        kernel_size = 13 #hvor stor kernel size er riktig?
        kernel = np.ones(kernel_size)/kernel_size
        y_s = np.convolve(y, kernel, mode='same')
        y_r[10:-10] = y_s[10:-10]
        return y_r

    def fit(x, y, expected, file, noise):
        params, cov = curve_fit(finder.bimodal_der, np.array([x[150:-150], noise[150:-150]]), y[150:-150], expected, bounds=((0, 0, 0, -np.inf), (np.inf, np.inf, np.inf, np.inf)))
        sigma = np.sqrt(np.diag(cov))
        d = finder.der(x, y)
        mid = int(np.average(np.where(d == np.nanmin(d))[0]))
        # plt.plot(x[150:-150], finder.bimodal_der(np.array([x, noise]), *params)[150:-150])
        plt.plot(x[150:-150], y[150:-150], "g")#NB: hvor stort område skal vi sammenlike, se på r2 kalling
        # plt.legend(["fit", "data"])
        # plt.plot(x[mid - 228], y[mid - 228], 'ro')
        # plt.plot(x[mid - 104], y[mid - 104], 'ro')
        # plt.plot(x[mid + 81], y[mid + 81], 'ro')
        # plt.plot(x[mid + 204], y[mid + 204], 'ro')
        plt.title(file)
        plt.grid()
        plt.xlabel('mT')
        plt.ylabel('Intensity')
        plt.show()
        # xyr = finder.xy_ratio(x, finder.bimodal_der(np.array([x, noise]), *params))
        xyr = finder.xy_ratio(x, y)



        high = np.amax(finder.bimodal_der(np.array([x, noise]), *params)[50:-50]) - np.amin(finder.bimodal_der(np.array([x, noise]), *params)[50:-50])

        return params, sigma, finder.r_sqr(finder.bimodal_der(np.array([x, noise]), *params)[75:-75], y[75:-75]), high, xyr #hvilket intervall skal vi se på?

    # def xy_ratio(x, y):
    #     for i in argrelextrema(y[500:650], np.less, order=50):
    #         print(float(x[500 + i]), float(y[500 + i]))
    #         min = float(y[500 + i])
    #     for i in argrelextrema(y[350:500], np.greater, order=10):
    #         print(np.average(x[350 + i]), np.average(y[350 + i]))
    #         y_r = np.average(y[350 + i]) - min
    #     for i in argrelextrema(y[225:325], np.greater, order=2):
    #         print(np.average(x[225 + i]), np.average(y[225 + i]))
    #         x_r = np.average(y[225 + i]) - min
    #     return x_r/y_r

    def xy_ratio(x, y):
        d = finder.der(x, y)
        mid = int(np.average(np.where(d == np.nanmin(d))[0]))


        x_1 = np.average(y[mid - 233: mid - 223])
        x_2 = np.average(y[mid + 199: mid + 209])
        y_1 = np.average(y[mid - 109: mid - 99 ])
        y_2 = np.average(y[mid + 76: mid + 86 ])
        y_r = y_1 - y_2
        x_r = x_1 - x_2

        return x_r/y_r



    def reader(file, expected, noise, k):
        f = open(file, "r")

        x = []
        y = []

        for line in f:
            x.append(float(line.split()[0]))
            y.append(float(line.split()[1]))

        x = np.array(x)
        y = np.array(y)

        if k == 1:
            params, sigma, r2, high, xyr = finder.fit(x, y, expected, file, noise)
            return params, sigma, r2, high, x, y, xyr
        else:
            return x, y

    def background(back):
        bac = []
        for file in back:
            b = open(file, "r")
            ba = []
            for line in b:
                ba.append(float(line.split()[1]))
            ba = np.array(ba)
            bac.append(ba)
        signal = np.zeros(len(ba))
        for i in range(len(bac)):
            signal += bac[i]


        return signal/len(bac)

    def conf_plot(x, y, type):
        p_x = np.arange(np.min(x), np.max(x) + 1, 1)

        if type==1:
            alpha, beta = np.polyfit(x, y, type)
            print(alpha, beta)
            p = 1
            p_y = alpha*x + beta
            print(finder.r_sqr(p_y, y))
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

            print(p_y)
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


        plt.xlabel('log(mp)')
        plt.ylabel('x/y ratio')

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
        print(finder.r_sqr(alpha*x + beta, y))
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

        plt.xlabel('log(mp)')
        plt.ylabel('x/y ratio')

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





db6 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\ivar_proton2\6\*")
db10 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\ivar_proton2\10\*")
db14 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\ivar_proton2\14\*")
db18 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\ivar_proton2\18\*")
dB6 = []
dB10 = []
dB14 = []
dB18 = []

for i in range(len(db6)):
    dB6.append(glob.glob(db6[i] + "\\*.dat"))
    dB10.append(glob.glob(db10[i] + "\\*.dat"))
    dB14.append(glob.glob(db14[i] + "\\*.dat"))
    dB18.append(glob.glob(db18[i] + "\\*.dat"))


dB6 = np.array(dB6, dtype='object')
dB10 = np.array(dB10, dtype='object')
dB14 = np.array(dB14, dtype='object')
dB18 = np.array(dB18, dtype='object')

for f in dB10[4]:
    x, y = finder.reader(f, 0, 0, k=0)
    print(f)
    plt.plot(x, y)
    plt.grid()
    plt.xlabel('mT')
    plt.show()
plt.plot(x, finder.background(dB10[0]))
plt.grid()
plt.xlabel('mT')
plt.show()
print(x[9] - x[1])

# trenger ikke den vide highw nå
expected = (20000, 3520, 4000, -1000)
# trenger ikke den vide highw nå
dB6 = np.delete(dB6, 2)
dB10 = np.delete(dB10, 2)
dB14 = np.delete(dB14, 2)
dB18 = np.delete(dB18, 2)

#kjører fit og skriver params til dokument
a1 = []
a2 = []
h = []
P = []
Y = []
columns = ["a1", "mu1", "a2", "bias"]
rows = ["s1-1", "s1-2", "s1-3", "s1-4", "s1-5", "s2-1", "s2-2", "s2-3", "s2-4", "s2-5", "s3-1", "s3-2", "s3-3", "s3-4", "s3-5", "s4-1", "s4-2", "s4-3", "s4-4", "s4-5", "s5-1", "s5-2", "s5-3", "s5-4", "s5-5"]


# doc = open("take2_RES18", "w")
unirad =finder.background(dB10[0])
print(unirad)

for f in np.concatenate(dB18[2:]):
    params, sigma, r2, high, x, y, xyr = finder.reader(f, expected, unirad, k=1)
    print(r2)
    print(params)
    a1.append(params[0])
    a2.append(params[2])
    h.append(high)
    P.append(params)
    Y.append(y)

# plt.plot(x, finder.background(dB10[0]))
# plt.show()

# plt.show()
# doc.write(str(pd.DataFrame(P, rows, columns)))
# doc.close()

pos = np.array([1, 2, 3, 4]*5)

plt.plot(pos, a1, "o")
plt.grid()
plt.title("a1")
plt.show()
plt.plot(pos, a2, "o")
plt.grid()
plt.title("a2")
plt.show()

data1 = np.zeros(len(x))
data3 = np.zeros(len(x))
peak1 = 0
peak3 = 0
for i in range(5):
    data1 = data1 + Y[i*4]
    peak1 = peak1 + h[i*4]
    data3 = data3 + Y[i*4 + 2]
    peak3 = peak3 + h[i*4 + 2]

data1 = data1/4
peak1 = peak1/4
data3 = data3/4
peak3 = peak3/4

plt.plot(x, data1/peak1)
plt.plot(x, data3/peak3)
plt.grid()
plt.legend(["1 (normalized with peak3/peak1)","3"], loc='upper right')
plt.title("With background")
plt.show()

plt.plot(x, (data1 - unirad)/(np.amax(data1 - unirad) - np.amin(data1 - unirad))) #2.0423253847667993
plt.plot(x, (data3 - unirad)/(np.amax(data3 - unirad) - np.amin(data3 - unirad)))
plt.grid()
plt.legend(["1 (normalized with peak3/peak1)","3"], loc='upper right')
plt.title("Without background")
plt.show()


# M = []
# T1 = []
# T2 = []
# B1 = []
# B2 = []
#
#
# dose_H = [dB6[1], dB10[1], dB14[1], dB18[1]]
# for f in np.concatenate(dose_H):
#     x, y = finder.reader(f,0, 0,k=0)
#     d = finder.der(x, y)
#     mid = int(np.where(d == np.nanmin(d))[0])
#     print(mid)
#     t1 = int(np.average(np.where(y[225:325] == np.amax(y[225:325]))[0]) + 225)
#     t2 = int(np.average(np.where(y[375:500] == np.amax(y[375:500]))[0]) + 375)
#     b1 = int(np.average(np.where(y[500:650] == np.amin(y[500:650]))[0]) + 500)
#     b2 = int(np.average(np.where(y[650:775] == np.amin(y[650:775]))[0]) + 650)
#
#
#     plt.plot(x, y)
#     plt.plot(x[mid], y[mid], 'go')
#     # plt.plot(x[t1], y[t1], 'ro') #-228
#     plt.plot(x[mid - 228], y[mid - 228], 'bo')
#     # plt.plot(x[t2], y[t2], 'ro') #-104
#     plt.plot(x[mid - 104], y[mid - 104], 'bo')
#     # plt.plot(x[b1], y[b1], 'ro') #81
#     plt.plot(x[mid + 81], y[mid + 81], 'bo')
#     # plt.plot(x[b2], y[b2], 'ro') #204
#     plt.plot(x[mid + 204], y[mid + 204], 'bo')
#     plt.plot(x, d)
#     plt.grid()
#     plt.title(f)
#     plt.show()
#     M.append(mid)
#     T1.append(t1)
#     T2.append(t2)
#     B1.append(b1)
#     B2.append(b2)
#
# M = np.array(M)[4:8]
# T1 = np.array(T1)[4:8]
# T2 = np.array(T2)[4:8]
# B1 = np.array(B1)[4:8]
# B2 = np.array(B2)[4:8]
#
# print(np.average(T1 - M))
# print(np.average(T2 - M))
# print(np.average(B1 - M))
# print(np.average(B2 - M))

#dette er avstanden fra vendepunkt til ekstremaler i 10 db kurver (de beste). Definerer dette som global avstand for disse punktene, ønsker videre å bruke dette til å finne xy ratio hvor vi gir noen punkter slingring på hver side og tar gjennomsnitt







#looking for dB patterns

dB = [np.concatenate(dB6[2:]), dB6[0], np.concatenate(dB10[2:]), dB10[0], np.concatenate(dB14[2:]), dB14[0], np.concatenate(dB18[2:]), dB18[0]]
# print(dB)
dba1 = []
dba2 = []
ratio = []
R2 = []

for i in range(4):
    a1 = []
    a2 = []
    xy = []
    l = dB[2*i]
    unirad = dB[2*i + 1]
    for f in l:
        # params, sigma, r2, high, x, y, xyr = finder.reader(f, expected, finder.background(unirad), k=1)
        # R2.append(r2)
        # print(r2)
        # a1.append(params[0])
        # a2.append(params[2])
        x, y = finder.reader(f, 0, 0, k=0)
        # print(len(x))
        xyr = finder.xy_ratio(x, y)
        xy.append(xyr)
    ratio.append(xy)
    # dba1.append(a1)
    # dba2.append(a2)
#
# f1_1 = []
# f1_2 = []
# f2_1 = []
# f2_2 = []
# f3_1 = []
# f3_2 = []
# f4_1 = []
# f4_2 = []

# print(np.average(R2))
R1 = []
R2 = []
R3 = []
R4 = []

for n in range(4):
    for i in range(5):
#         # f1_1.append(dba1[n][i*4])
#         # f1_2.append(dba2[n][i*4])
#         # f2_1.append(dba1[n][i*4 + 1])
#         # f2_2.append(dba2[n][i*4 + 1])
#         # f3_1.append(dba1[n][i*4 + 2])
#         # f3_2.append(dba2[n][i*4 + 2])
#         # f4_1.append(dba1[n][i*4 + 3])
#         # f4_2.append(dba2[n][i*4 + 3])
#

        R1.append(ratio[n][i*4])
        R2.append(ratio[n][i*4 + 1])
        R3.append(ratio[n][i*4 + 2])
        R4.append(ratio[n][i*4 + 3])
#
#
#
dbx = np.array([6]*5 + [10]*5 + [14]*5 + [18]*5)
mp = np.array([50.45]*5 + [20.08]*5 + [7.996]*5 + [3.182]*5) #tall hentet fra .par filer, men litt usikekr på om dette faktisk er micorwave power
#
# # plt.plot(mp, f1_1, 'o')
# # plt.plot(mp, f2_1, 'o')
# # plt.plot(mp, f3_1, 'o')
# # plt.plot(mp, f4_1, 'o')
# # plt.xlabel('mp')
# # plt.ylabel('a1')
# # plt.legend(['f1', 'f2', 'f3', 'f4'])
# # plt.grid()
# # plt.show()
# #
# # plt.plot(mp, f1_2, 'o')
# # plt.plot(mp, f2_2, 'o')
# # plt.plot(mp, f3_2, 'o')
# # plt.plot(mp, f4_2, 'o')
# # plt.xlabel('mp')
# # plt.ylabel('a2')
# # plt.legend(['f1', 'f2', 'f3', 'f4'])
# # plt.grid()
# # plt.show()

plt.plot(np.log(mp), R1, 'o')
plt.plot(np.log(mp), R2, 'o')
plt.plot(np.log(mp), R3, 'o')
plt.plot(np.log(mp), R4, 'o')
plt.legend(['film 1', 'film 2', 'film 3', 'film 4'])
plt.xlabel('log(mp)')
plt.ylabel('x/y ratio')
# plt.xlim([3, 100])
# plt.xscale("log")
plt.grid()
plt.show()

print('Film 1')
finder.conf_plot2(np.log(mp), R1)
print('Film 2')
finder.conf_plot2(np.log(mp), R2)
print('Film 3')
finder.conf_plot2(np.log(mp), R3)
print('Film 4')
finder.conf_plot2(np.log(mp), R4)





# finder.conf_plot(np.log(mp), R1, 1)
# finder.conf_plot(np.log(mp), R2, 1)
# finder.conf_plot(np.log(mp), R3, 1)

#
#
# x1, y1 = finder.reader(dB14[1][0], 0, 0, k=0)
# x2, y2 = finder.reader(dB14[1][1], 0, 0, k=0)
# x3, y3 = finder.reader(dB14[1][2], 0, 0, k=0)
# x4, y4 = finder.reader(dB14[1][3], 0, 0, k=0)
#
#
# plt.plot(x1, y1/(np.amax(y1) - np.amin(y1)))
# plt.plot(x2, y2/(np.amax(y2) - np.amin(y2)))
# plt.plot(x3, y3/(np.amax(y3) - np.amin(y3)))
# plt.plot(x4, y4/(np.amax(y4) - np.amin(y4)))
# plt.legend(['film 1', 'film 2', 'film 3', 'film 4'])
# plt.xlabel('mT')
# plt.grid()
# plt.show()
#
# #lager bragg kurve plott
#
# a1 = []
# a2 = []
# h = []
# P = []
# Y = []
# unirad =finder.background(dB10[0])
#
# for f in np.concatenate(dB10[2:7]):
#     params, sigma, r2, high, x, y, xyr = finder.reader(f, expected, unirad, k=1)
#     print(r2)
#     # print(params)
#     a1.append(params[0])
#     a2.append(params[2])
#     h.append(high)
#     P.append(params)
#     Y.append(y)
#
# plt.plot(x, unirad)
# plt.grid()
# plt.xlabel('mT')
# plt.show()
#
# pos = np.array([1, 2, 3, 4]*5)
#
# plt.plot(pos, a1, "o")
# plt.grid()
# plt.title("a1")
# plt.xlabel('Film nr.')
# plt.ylabel('amp')
# plt.show()
# plt.plot(pos, a2, "o")
# plt.grid()
# plt.title("a2")
# plt.xlabel('Film nr.')
# plt.ylabel('amp')
# plt.show()
