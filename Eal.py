import numpy as np
import matplotlib.pyplot as plt
import glob





def reader(file):
    f = open(file, "r")
    x = []
    y = []

    for line in f:
        x.append(float(line.split()[0]))
        y.append(float(line.split()[1]))

    x = np.array(x)
    y = np.array(y)

    n = np.where(y[120:260]==np.nanmax(y[120:260]))
    n = int(np.average(n[0]))
    m = np.where(y[260:400]==np.nanmin(y[260:400]))
    m = int(np.average(m[0]))

    plt.plot(x, y)
    plt.plot(x[120 + n], y[120 + n], 'ro')
    plt.plot(x[260 + m], y[260 + m], 'ro')
    plt.grid()
    plt.xlabel('mT')
    plt.show()

    amp = y[120 + n] - y[260 + m]
    return amp


#
# mar06 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\Xrays_Ala_Cu_06 Mar 23\*.dat")
#
#
# VAL = []
# Lcon = 0
# ENL = 0
# new = []
# for f in mar06:
#     print(f)
#     F = f.split()
#     con = F[2]
#     REST = F[-1]
#     R = REST.split('_')
#     W = float(R[3])
#     # print(W)
#     EN = R[2]
#     print(EN)
#     if con == Lcon and EN == ENL:
#         print('equal consentration')
#         print('equal energy')
#         peaks = reader(f)/W
#         new.append(peaks)
#         Lcon = con
#         ENL = EN
#     else:
#         # print(f)
#         VAL.append(new)
#         print('new')
#         new = []
#         peaks = reader(f)/W
#         new.append(peaks)
#         Lcon = con
#         ENL = EN
#
# VAL.append(new)
# VAL.pop(0)

mar07 = glob.glob(r"C:\Users\Hp\Desktop\skole\msc\Xrays_Ala_Cu_07 Mar 23\*.dat")

VAL2 = []
Lcon2 = 0
ENL2 = 0
new2 = []
for f in mar07:
    print(f)
    F = f.split()
    con = F[2]
    REST = F[-1]
    R = REST.split('_')
    W = float(R[3])
    # print(W)
    EN = R[2]
    print(EN)
    if con == Lcon2 and EN == ENL2:
        # print('equal consentration')
        # print('equal energy')
        peaks = reader(f)/W
        new2.append(peaks)
        Lcon2 = con
        ENL2 = EN
    else:
        # print(f)
        VAL2.append(new2)
        print('new')
        new2 = []
        peaks = reader(f)/W
        new2.append(peaks)
        Lcon2 = con
        ENL2 = EN

VAL2.append(new2)
VAL2.pop(0)


AVG = np.zeros(len(VAL))
for i in range(len(VAL)):
    AVG[i] = np.average(VAL[i])
AVG2 = np.zeros(len(VAL2))
for i in range(len(VAL2)):
    AVG2[i] = np.average(VAL2[i])


plt.plot(120, AVG[12], color='tab:orange', marker='o', label='Pure Alanine')
plt.plot(96.6, AVG[13], color='tab:orange', marker='o')
plt.plot(72.2, AVG2[-1], color='tab:orange', marker='o')
plt.plot(48.8, AVG2[-2], color='tab:orange', marker='o')
plt.plot(25.09, AVG2[-3], color='tab:orange', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[-3], AVG2[-2], AVG2[-1], AVG[13], AVG[12]], color='tab:orange')

plt.plot(120, AVG[0], color='tab:green', marker='o',label='0.5 mol%')
plt.plot(96.6, AVG[1], color='tab:green', marker='o')
plt.plot(72.2,AVG2[2] , color='tab:green', marker='o')
plt.plot(48.8, AVG2[1] , color='tab:green', marker='o')
plt.plot(25.09,AVG2[0], color='tab:green', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[0], AVG2[1], AVG2[2], AVG[1], AVG[0]], 'g')

plt.plot(120, AVG[4], color='tab:red', marker='o', label='1.0 mol%')
plt.plot(96.6, AVG[5], color='tab:red', marker='o')
plt.plot(72.2, AVG2[8], color='tab:red', marker='o')
plt.plot(48.8, AVG2[7], color='tab:red', marker='o')
plt.plot(25.09,AVG2[6] , color='tab:red', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[6], AVG2[7], AVG2[8], AVG[5], AVG[4]], color='tab:red')

plt.plot(120, AVG[2], color='tab:purple', marker='o', label='1.5 mol%')
plt.plot(96.6, AVG[3], color='tab:purple', marker='o')
plt.plot(72.2,AVG2[5] , color='tab:purple', marker='o')
plt.plot(48.8, AVG2[4], color='tab:purple', marker='o')
plt.plot(25.09, AVG2[3], color='tab:purple', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[3], AVG2[4], AVG2[5], AVG[3], AVG[2]], color='tab:purple')

plt.plot(120, AVG[8], color='tab:brown', marker='o', label='2.0 mol%')
plt.plot(96.6, AVG[9], color='tab:brown', marker='o')
plt.plot(72.2, AVG2[14], color='tab:brown', marker='o')
plt.plot(48.8, AVG2[13], color='tab:brown', marker='o')
plt.plot(25.09, AVG2[12], color='tab:brown', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[12], AVG2[13], AVG2[14], AVG[9], AVG[8]], color='tab:brown')

plt.plot(120, AVG[6], color='tab:pink', marker='o', label='2.5 mol%')
plt.plot(96.6, AVG[7], color='tab:pink', marker='o')
plt.plot(72.2, AVG2[11], color='tab:pink', marker='o')
plt.plot(48.8, AVG2[10], color='tab:pink', marker='o')
plt.plot(25.09,AVG2[9] , color='tab:pink', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[9], AVG2[10], AVG2[11], AVG[7], AVG[6]], color='tab:pink')

plt.plot(120, AVG[10], color='tab:gray', marker='o', label='3.0 mol%')
plt.plot(96.6, AVG[11], color='tab:gray', marker='o')
plt.plot(72.2,AVG2[17] , color='tab:gray', marker='o')
plt.plot(48.8, AVG2[16], color='tab:gray', marker='o')
plt.plot(25.09, AVG2[15], color='tab:gray', marker='o')
plt.plot([25.09, 48.8, 72.2, 96.6, 120], [AVG2[15], AVG2[16], AVG2[17], AVG[11], AVG[10]], color='tab:gray')


plt.legend(bbox_to_anchor=(1.0, 1.0))
leg = plt.gca().get_legend()
ltext = leg.get_texts()
plt.setp(ltext, fontsize=10)
plt.grid()
plt.xlabel('keV')
plt.ylabel('Relative amplitude (corrected for weight)')
plt.show()
