import numpy as np
import matplotlib.pyplot as plt
import path
import csv

def read_aero_coefficients_from_csv(filename, N):
    """
    Чтение аэродинамических коэффициентов из CSV файла
    """
    # Инициализация массивов нулями
    Cwv = np.zeros(N)
    Cww = np.zeros(N)
    Cwb = np.zeros(N)
    Cvv = np.zeros(N)
    Cvb = np.zeros(N)
    Cvw = np.zeros(N)
    Wind = np.zeros(N)
    
    try:
        with open(filename, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            
            i = 0
            for row in reader:
                if i >= N:
                    break
                    
                Cvw[i] = float(row['Cyw'])
                Cww[i] = float(row['Cww'])
                Cwb[i] = float(row['Cwb'])
                Cvv[i] = float(row['Cyy'])
                Cvb[i] = -float(row['Cbs'])
                Cvw[i] = float(row['Cyw'])
                Wind[i] = float(row['wind']) 
                
                i += 1
                
        print(f"Успешно прочитано {i} строк из {filename}")
        
    except FileNotFoundError:
        print(f"Файл {filename} не найден")
    except Exception as e:
        print(f"Ошибка при чтении файла {filename}: {e}")
    
    return Cwv, Cww, Cwb, Cvv, Cvb, Cvw, Wind

# Основной код
N = 3840

# Инициализация массивов
Cwv, Cww, Cwb, Cvv, Cvb, Cvw, Wind = read_aero_coefficients_from_csv("output/dynamic_coefs.csv", N)

# Остальные массивы (не из файла)
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

mass_st = 14400
l_1 = 23720

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

    if t<100:
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