import numpy as np
import matplotlib.pyplot as plt
import path


N = 3840
Cwv = np.empty(N)
Cww = np.empty(N)
Cwb = np.empty(N)
Cvv = np.empty(N)
Cvb = np.empty(N)
Cvw = np.empty(N)
Wind = np.empty(N)
Ms = np.empty(N)
W1 = np.empty(N)
W2 = np.empty(N)
W3 = np.empty(N)
W4 = np.empty(N)
W5 = np.empty(N)
Cy1 = np.empty(N)
Cy2 = np.empty(N)
Cy3 = np.empty(N)
Cy4 = np.empty(N)
Cy5 = np.empty(N)
Cws = np.empty(N)
Cw1 = np.empty(N)
Cw2 = np.empty(N)
Cw3 = np.empty(N)
Cw4 = np.empty(N)
Cw5 = np.empty(N)
Cws = np.empty(N)

mass_st = 14400
l_1 = 23720



with open(path.root_path + "Aero.txt", "r") as Aero:
    for i in range(N):

        try:
            # Читаем и обрабатываем каждую строку
            Cwv[i] = float(Aero.readline().strip())
            Cww[i] = float(Aero.readline().strip())
            Cwb[i] = float(Aero.readline().strip())
            Cvv[i] = float(Aero.readline().strip())
            Cvb[i] = -float(Aero.readline().strip())
            Cvw[i] = float(Aero.readline().strip())
            Wind[i] = float(Aero.readline().strip())
        except ValueError as e:
            # Выводим ошибку и строку, которая вызвала её
            #print(f"Ошибка преобразования строки на итерации {i}: {e}")
            #print(f"Непрочитанное значение: '{Aero.readline().strip()}'")  # Выводим непрочитанное значение
            # Задаем значение по умолчанию
            Cwv[i] = 0.0
            Cww[i] = 0.0
            Cwb[i] = 0.0
            Cvv[i] = 0.0
            Cvb[i] = 0.0
            Cvw[i] = 0.0
            Wind[i] = 0.0

h = 0.1
uc = 0
duc = 0
dduc = 0
v = 0
y = 0
dv = 0
w = 0
dw = 0
ddw = 0
t = 0
t1 = 0.38
t2 = 0.04

X = []
Y = []
Z = []
K = []
A = []
B = []
V = []
W = []

zerX = []
zerY = []

a0 = 1.32
a1 = 1.7
a2 = 0.0004
a3 = 10*a2


c = 0

for i in range(N):

    dv =  Cvv[i] * Wind[i] - Cvw[i] * w - Cvv[i] * v - Cvb[i] * uc
    ddw = Cwv[i] * Wind[i] - Cww[i] * w - Cwv[i] * v - Cwb[i] * uc

    dduc = (-t1*duc - uc + a0*w + a1*dw + a2*y * a3*v)/t2
    #dduc = signal(dw, Kpw, Kiw, Kdw)

    if (abs(t-0.2)<0.0001) or (abs(t-30.0)<0.0001) or (abs(t-50.0)<0.0001):
            print("y={:.3e} ({})".format(y, t))
            print("v={:.3e} ({})".format(v, t))
            print("a={:.3e} ({})".format(dv, t))
            print("w={:.3e} ({})".format(w, t))
            print("dw={:.3e} ({})".format(dw, t))
            print("ddw={:.3e} ({})".format(ddw, t))
            print("b={:.3e} ({})".format(uc, t))


            print(f"Iteration {i}:")
            print(f"  Cvv[i]: {Cvv[i]:.3e}")
            print(f"  Wind[i]: {Wind[i]:.3e}")
            print(f"  Cvw[i]: {Cvw[i]:.3e}")
            print(f"  w: {w:.3e}")
            print(f"  v: {v:.3e}")
            print(f"  Cvb[i]: {Cvb[i]:.3e}")
            print(f"  uc: {uc:.3e}")
            print(f"  Cwv[i]: {Cwv[i]:.3e}")
            print(f"  Cww[i]: {Cww[i]:.3e}")
            print(f"  Cwb[i]: {Cwb[i]:.3e}")
            print(f"  dv: {dv:.3e}")
            print(f"  ddw: {ddw:.3e}")
            print(dv-(Cvv[i] * Wind[i] - Cvw[i] * w - Cvv[i] * v - Cvb[i] * uc))
            print()
    v += h * dv
    y += h * v

    dw += h * ddw
    w += h * dw
    t += h
    duc += h*dduc
    uc += h*duc

    if uc > 7/57.3:
        uc = 7/57.3

    if uc < -7/57.3:
        uc = -7/57.3

    if uc > 0 and uc < 0.3/57.3:
       uc = 0.3/57.3

    if uc < 0 and uc > -0.3/57.3:
       uc = -0.3/57.3

    if t<52:
        Y.append(float(v))
        X.append(float(t))

        Z.append(float(t))
        K.append(float(w*57.3))

        A.append(float(t))
        B.append(float(y))

        zerY.append(float(0))
        zerX.append(float(t))

        V.append(float(uc*57.3))
        W.append(float(t))

    # if c==0:
    #     if (t >= 365.2):
        #if (abs(Cwb[i]-0.0)>0.00001):

        #     c+=1
    #print(dduc*57.3, "____ dduc|duc ____", duc*57.3)

plt.subplot(4, 1, 1)
plt.plot(X, Y)
plt.ylabel('Скорость(t), м/c',color='gray')
plt.plot(zerX, zerY)
plt.grid(True)

plt.subplot(4, 1, 2)
plt.plot(Z, K)
plt.ylabel('Угол (t), град',color='gray')
plt.plot(zerX, zerY)
plt.grid(True)

plt.subplot(4, 1, 3)
plt.plot(A, B)
plt.ylabel('Перемещение(t), м',color='gray')
plt.plot(zerX, zerY)
plt.grid(True)

plt.subplot(4, 1, 4)
plt.plot(W, V)
plt.ylabel('Поворот ОУ(t), град',color='gray')
plt.plot(zerX, zerY)
plt.grid(True)

plt.show()