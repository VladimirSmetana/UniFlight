import csv
import matplotlib.pyplot as plt
import math as m
import path
import rocket_parser as rp
import constants

# Инициализация списков
numeric = []
length = []
read_mass = []
mass = []
stiffness = []

L = [4.73,  7.853, 10.996, 14.137, 17.279]

def absmax(iterable):
    return max(iterable, key=abs)

def calculate_sum(base):
    sum = [0] * len(base)
    for i in range(len(base)):
        if i == 0:
            sum[i] = base[i]
        else:
            sum[i] = sum[i-1] + base[i - 1] + base[i]
    return sum

def calculate_multi(one, second):
    return [a * b for a, b in zip(one, second)]


# Чтение данных из CSV-файла
n = 0
length_max = 65
with open(path.root_path + 'rocket_body.csv', newline='') as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')  # Используем запятую как разделитель
    next(spamreader)  # Пропускаем заголовок
    for row in spamreader:
        numeric.append(n)
        length.append(float(row[0]))  # Преобразуем значения в float и добавляем в список
        read_mass.append(float(row[1]))
        stiffness.append(float(row[2]))
        n+=1
        if n > length_max:
            break


parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")
step = parser.get_interstep()
work_time = parser.get_work_time()



def changed_mass(current_time):
    block_number = parser.get_block_number()

    delta_level_fu = parser.get_delta_level_fu()
    delta_mass_fu = parser.get_delta_mass_fu()
    sector_range_fu = parser.get_sector_index_fu()

    delta_level_ox = parser.get_delta_level_ox()
    delta_mass_ox = parser.get_delta_mass_ox()
    sector_range_ox = parser.get_sector_index_ox()
    
    mass_t = read_mass.copy()

    for j in range(block_number):
        # Обработка сектора окислителя
        sector_number_ox = sector_range_ox[j][1] - sector_range_ox[j][0]
        for k in range(*sector_range_ox[j]):
            time = 0
            while time < work_time[j] / sector_number_ox:
                if time >= current_time:
                    break
                mass_t[k] -= (delta_mass_ox[j] / 1000 * step)
                mass_t[k] = max(mass_t[k], 0)  # Не допускаем отрицательной массы
                time += step

        # Обработка сектора топлива
        sector_number_fu = sector_range_fu[j][1] - sector_range_fu[j][0]
        for k in range(*sector_range_fu[j]):
            time = 0
            while time < work_time[j] / sector_number_fu:
                if time >= current_time:
                    break
                mass_t[k] -= (delta_mass_fu[j] / 1000 * step)
                mass_t[k] = max(mass_t[k], 0)  # Не допускаем отрицательной массы
                time += step

    return mass_t


# Начальные параметры это учитывают
ti = 0.0
ver_mass_vector = []
time_vector = []
frec_vector = []
while ti < work_time[0]:
    ver_mass_vector.append(changed_mass(ti))
    time_vector.append(ti)
    ti += step

start_color = [0.68, 0.85, 0.9]
end_color = [0, 0, 0.55]   
total_iterations = len(ver_mass_vector)

def interpolate_color(start_color, end_color, i, total):
    return [
        start_color[j] + (end_color[j] - start_color[j]) * i / (total - 1)
        for j in range(3)
    ]

en = 0
for mass in ver_mass_vector:
    def delta_vector(previous, actual):
        f_12 = [x ** 2 for x in previous]
        mf_12 = calculate_multi(f_12, mass)
        sum_mf_12 = calculate_sum(mf_12)
        f1_f20 = calculate_multi(previous, actual)
        f1_f20_mass = calculate_multi(f1_f20, mass)
        f1_f20_mass_summ = calculate_sum(f1_f20_mass)
        delta12 = - f1_f20_mass_summ[-1]/sum_mf_12[-1]
        delta12f = [x * delta12 for x in previous]
        return delta12f

    m_N = [a * b for a, b in zip(numeric, mass)]

    sum_m = calculate_sum(mass)
    sum_m_N = calculate_sum(m_N)
    Nmid = sum_m_N[-1]/sum_m[-1]

    N_Nm  = [x - Nmid for x in numeric]
    N_Nm2 = [x ** 2 for x in N_Nm]

    preIn = calculate_multi(N_Nm2, mass)
    In = calculate_sum(preIn)

    rocket_length = 0
    rocket_mass = 0
    for i in range(len(length)):
        rocket_length += length[i]
        rocket_mass   += mass[i]

    a = [x/rocket_length for x in L]
    chJ_cosJ = [m.cosh(x)-m.cos(x) for x in L]
    shJ_sinJ = [m.sinh(x)-m.sin(x) for x in L]
    sinJ_shJ = [m.sin(x)-m.sinh(x) for x in L]
    Y = [a / b for a, b in zip(chJ_cosJ, sinJ_shJ)]


    f_zero = [0] * 5
    f_stiffness = [0] * 5
    f_mass = [0] * 5
    w_calc = [0] * 5
    fi     = [0] * 5
    for i in range(len(f_zero)):
        f_zero[i] = list(((m.sin(a[i]*x)+m.sinh(a[i]*x))*Y[i]+(m.cos(a[i]*x)+m.cosh(a[i]*x)))/2 for x in numeric)

    w_zero = [m.sqrt(max(stiffness)/(rocket_mass*(10**3)/rocket_length*pow(rocket_length,4)))*(x**2)/(2*m.pi) for x in L]

    def calculate_form(index):

        f_start = f_zero[index]
        tolerance = 1e-8
        while (True):
            m_f1 = calculate_multi(mass, f_start)
            sum_m_f1 = calculate_sum(m_f1)

            value_6_11 = calculate_multi(m_f1, N_Nm)
            sum_value_6_11 = calculate_sum(value_6_11)

            D1 = - sum_value_6_11[-1]/In[-1]
            D2 = - sum_m_f1[-1]/sum_m[-1]

            D1_6 = [x*D1 for x in N_Nm]
            D2_15 = [D2+x  for x in D1_6]
            accumulated_delta = [0] * len(f_start)
            if index > 0:
                for i in range(index):
                    newin = delta_vector(f_mass[index - 1 - i], f_start)
                    accumulated_delta = [a  + b  for a, b in zip(accumulated_delta, newin)]

            f1_16 = [a + b + c for a, b, c in zip(D2_15, f_start, accumulated_delta)]
            f_mass[index] = [x/max(f1_16) for x in f1_16]

            m_f1 = calculate_multi(mass, f_mass[index])
            sum_m_f1 = calculate_sum(m_f1)
            double_sum_m_f1 = calculate_sum(sum_m_f1)
            dm1 = [-x*double_sum_m_f1[-1]/numeric[-1] for x in numeric]
            M1x = [a + b for a, b in zip(dm1, double_sum_m_f1)]
            M1x_E = [a/b if b != 0 else 0 for a, b in zip(M1x, stiffness)]
            sum_M1x_E = calculate_sum(M1x_E)
            fi[index] = calculate_sum(sum_M1x_E)
            double_sum_M1x_E_mass = calculate_multi(fi[index], mass)
            summ_13 = calculate_sum(double_sum_M1x_E_mass)

            value_13_15 = [a * b for a, b in zip(double_sum_M1x_E_mass, N_Nm)]
            sum_13_15 = calculate_sum(value_13_15)

            D1 = - sum_13_15[-1]/In[-1]
            D2 = - summ_13[-1]/sum_m[-1]

            D1_15 = [x*D1 for x in N_Nm]

            D2_11 = [a + b + D2 for a, b in zip(fi[index], D1_15)]
            accumulated_delta = [0] * len(fi[index])
            if index > 0:
                for i in range(index):
                    newin = delta_vector(f_stiffness[index - 1 - i], fi[index])
                    accumulated_delta = [(a + b) for a, b in zip(accumulated_delta, newin)]

            D2_11 = [a + b  for a, b in zip(D2_11, accumulated_delta)]

            f_stiffness_res = [x/absmax(D2_11) for x in D2_11]


            m_f1 = calculate_multi(mass, f_stiffness_res)
            sum_m_f1 = calculate_sum(m_f1)
            double_sum_m_f1 = calculate_sum(sum_m_f1)
            dm1 = [-x*double_sum_m_f1[-1]/numeric[-1] for x in numeric]
            M1x = [a + b for a, b in zip(dm1, double_sum_m_f1)]
            M1x2 = [x ** 2 for x in M1x]
            M1x2_E = [a/b if b != 0 else 0 for a, b in zip(M1x2, stiffness)]
            sum_M1x2_E = calculate_sum(M1x2_E)
            f_12 = [x ** 2 for x in f_stiffness_res]
            mf_12 = calculate_multi(f_12, mass)
            sum_mf_12 = calculate_sum(mf_12)
            w_calc[index] = m.sqrt(sum_mf_12[-1]/(sum_M1x2_E[-1]*1000.0*pow(length[-1]/2,4)))/(2*m.pi)

            f_start = f_stiffness_res
            if max(abs(f_stiffness_res[i] - f_start[i]) for i in range(len(f_start))) < tolerance:
                break

        return f_stiffness_res
    #################################################################



    w_femap = [11.86, 32.51, 60.66, 86.75, 124.37]
    color_pairs = [
        ([0.68, 0.85, 0.9], [0, 0, 1]),
        ([1, 0.68, 0.68], [1, 0, 0]),
        ([0.8, 0.8, 0.8], [0, 0, 0])
    ]

    for i in range(0, 3):
        f_stiffness[i] = calculate_form(i)
        # print("w["+str(i)+"] = " + str((w_calc[i])) + " / " + str((w_femap[i])) + " -> " + str(abs(m.floor((w_calc[i] - w_femap[i]) * 100 /w_femap[i]))) +" %")
        plt.plot(numeric, f_stiffness[i], color = interpolate_color(color_pairs[i][0], color_pairs[i][1], en, total_iterations))
                #  , label = [f'{i+1} Тон - {round(w_calc[i], 2)} Hz'])
    frec_vector.append(w_calc[0])

    plt.title('Расчет форм и частот колебаний', fontsize=16)
    plt.xlabel('Длина РН, м', fontsize=14)
    plt.ylabel('Форма', fontsize=14)
    plt.grid(True)
    plt.tight_layout()
    en+=1
    print(mass)

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Создаем кастомные элементы для легенды
custom_lines = [
    Line2D([0], [0], color='blue', lw=2),
    Line2D([0], [0], color='black', lw=2),
    Line2D([0], [0], color='red', lw=2)
]

plt.legend(custom_lines, ['1 Тон', '2 Тон', '3 Тон'])

plt.show()

plt.plot(time_vector, frec_vector)
plt.show()
