import rocket_parser as rp
import constants
import csv
import os
import path
import atmosphere as atmo
import aerodynamics as aero
import math
from scipy.integrate import solve_ivp
import attack

parser = rp.rocket_parser(path.rocket_lib + "falcon9.json")

alpha = attack.alpha(parser.attack_coefs[0],
                     parser.attack_coefs[1],
                     parser.work_time[0],
                     False)

def get_attack(vel, time):
    return alpha.calculate_alpha(vel, time)

def write_arrays_to_csv(filename, **arrays):
    if not arrays:
        raise ValueError("Array is required.")
    
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    
    headers = list(arrays.keys())
    max_length = min(len(arr) for arr in arrays.values())
    
    with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(headers)
        for i in range(max_length):
            row = [arrays[name][i] for name in headers]
            writer.writerow(row)
    
    print(f"Data was moved to '{filename}'.")


Cbs_list = []
Cyw_list = []
Cww_list = []
Cyy_list = []
Cwy_list = []
Cwb_list = []
Csb_list = []

attack_list = []
time_list = []
wind_list = []

class ballistics:
    def __init__(self, N, Y, vel, alt):
        self.N = N
        self.Y = Y
        self.vel = vel
        self.alt = alt
        
        self.G = aero.UnionStream()
        self.G.set_elnumber(parser.get_block_number()+1)
        self.G.set_diameter(parser.get_diameters())
        self.G.set_length(parser.get_part_length())

        self.thrust = 0
        self.mass = 0
        self.inertia = 0
        self.attack = 0
        self.dencity = 0
        self.dypressure = 0
        self.first_point = 0
        
        self.atm = atmo.atmosphere(self.alt)
        self.last_time = None

    def update_params(self, time):
        if self.last_time != time:
            self.thrust  =  parser.get_thrust_from_time(time)
            self.mass    =  parser.get_mass_from_time(time)
            self.inertia =  parser.get_inertia_from_time(time)
            self.center  =  parser.get_center_from_time(time) 
            
            # Проверка на None значения
            if self.center is None:
                print(f"Warning: center is None at time {time}, using default value 0.0")
                self.center = 0.0
            if self.thrust is None:
                print(f"Warning: thrust is None at time {time}, using default value 0.0")
                self.thrust = 0.0
            if self.mass is None:
                print(f"Warning: mass is None at time {time}, using default value 0.1")
                self.mass = 0.1
            if self.inertia is None:
                print(f"Warning: inertia is None at time {time}, using default value 0.1")
                self.inertia = 0.1
                
            self.attack  = get_attack(self.vel, time)*math.pi/180
            self.G.calculate_CXY(self.vel, self.alt, self.attack)

            self.atm = atmo.atmosphere(self.alt)
            self.dencity = self.atm.get_density()
            self.wind = self.atm.get_wind()
            
            # Проверка focus_position
            if not hasattr(self.G, 'focus_position') or self.G.focus_position is None:
                print(f"Warning: focus_position is None at time {time}, using default value 0.0")
                focus_pos = 0.0
            else:
                focus_pos = self.G.focus_position
                
            self.first_point  = abs(focus_pos - self.center)
            self.second_point = abs(parser.rocket_length - self.center)

            if self.alt > 90000:
                self.G.CX = 0
                self.G.CY = 0
                self.dencity = 0
            
            if self.alt < 0:
                self.alt = 0

            self.last_time = time
            self.dypressure = self.dencity * self.vel**2/2

            #output
            attack_list.append(self.attack*180/math.pi)
            time_list.append(time)
            wind_list.append(self.wind)
            Cbs_list.append(self.thrust*parser.thrust_ratio/self.mass) # как по рысканию
            Cyw_list.append(-(self.thrust+self.G.CY*self.dypressure*parser.maximum_area)/self.mass) # как по рысканию, но со знаком "-"
            Cww_list.append((-self.G.CY*self.dypressure*parser.maximum_area*self.first_point)/self.inertia) # как по рысканию, но со знаком "-"
            Cyy_list.append((self.G.CY*self.dypressure*parser.maximum_area)/(self.mass*self.vel)) # как по рысканию
            Cwy_list.append((self.G.CY*self.dypressure*parser.maximum_area*self.first_point)/self.inertia/self.vel) # как по рысканию
            Cwb_list.append(self.thrust*parser.thrust_ratio*self.second_point/self.inertia) # как по рысканию
            Csb_list.append(self.thrust*parser.thrust_ratio/self.inertia) # как по рысканию

    def delta_velocity(self, time):
        self.update_params(time)
        F_P = self.thrust * math.cos(self.attack)
        F_X = self.G.CX*self.dypressure*parser.maximum_area
        return (F_P - F_X)/self.mass - self.atm.get_AOG() * math.sin(self.Y)

    def delta_trajangle(self, time):
        self.update_params(time)
        F_P = self.thrust * math.sin(self.attack)
        F_Y = self.G.CY*self.dypressure*parser.maximum_area
        F_G = self.atm.get_AOG() * math.cos(self.Y) * (1 - self.vel**2 / (self.atm.get_AOG() * (constants.earth_radius + self.alt))) 
        return (F_P + F_Y)/(self.mass*self.vel) - F_G/self.vel
    
    def delta_polar(self, time):
        self.update_params(time)
        return (self.vel/(constants.earth_radius + self.alt))*math.cos(self.Y)
    
    def delta_altitude(self, time):
        self.update_params(time)
        return self.vel * math.sin(self.Y)

    def delta_longitude(self, time):
        self.update_params(time)
        return self.vel * math.cos(self.Y)
    
def output():
    rocketname = parser.name
    write_arrays_to_csv("output/"+rocketname+"_dynamic_coefs.csv",
                        time=time_list,
                        wind=wind_list,
                        Cbs=Cbs_list,
                        Cyw=Cyw_list,
                        Cww=Cww_list,
                        Cyy=Cyy_list,
                        Cwy=Cwy_list,
                        Cwb=Cwb_list,
                        Csb=Csb_list)

def system(t, vars):
    n, y, v, h, l = vars
    b = ballistics(n, y, v, h)
    return [
        b.delta_polar(t),
        b.delta_trajangle(t),
        b.delta_velocity(t),
        b.delta_altitude(t),
        b.delta_longitude(t)
    ]

def check_parser_data():
    """Проверяет данные парсера перед запуском расчета"""
    print("=== ПРОВЕРКА ДАННЫХ ПАРСЕРА ===")
    
    ft = parser.get_full_time()
    print(f"Полное время полета: {ft} с")
    
    # Проверяем ключевые параметры на разных временных точках
    test_times = [0, ft/2, ft-1]
    
    for t in test_times:
        print(f"\nВремя {t} с:")
        thrust = parser.get_thrust_from_time(t)
        mass = parser.get_mass_from_time(t)
        inertia = parser.get_inertia_from_time(t)
        center = parser.get_center_from_time(t)
        
        print(f"  Тяга: {thrust}")
        print(f"  Масса: {mass}")
        print(f"  Инерция: {inertia}")
        print(f"  Центр: {center}")
        
        if any(x is None for x in [thrust, mass, inertia, center]):
            print(f"  ⚠️  ВНИМАНИЕ: Найдены None значения!")
    
    print(f"\nДлина ракеты: {parser.rocket_length}")
    print(f"Максимальная площадь: {parser.maximum_area}")
    print(f"Коэффициент тяги: {parser.thrust_ratio}")

def main():
    """Основная функция с проверками"""
    print("ЗАПУСК БАЛЛИСТИЧЕСКОГО РАСЧЕТА")
    print("=" * 50)
    
    # Проверка данных
    check_parser_data()
    
    print("\n" + "=" * 50)
    print("ЗАПУСК ИНТЕГРИРОВАНИЯ...")
    
    try:
        ft = parser.get_full_time()
        h = parser.interstep
        t_span = (0, ft-1)
        y0 = [0, math.pi/2, 0.1, 0.1, 0.1]
        
        print(f"Временной диапазон: {t_span}")
        print(f"Шаг интегрирования: {h}")
        print(f"Начальные условия: {y0}")
        
        # Запуск интегрирования
        sol = solve_ivp(system, t_span, y0, method='RK45', max_step=h)
        
        print("\n=== РЕЗУЛЬТАТЫ ===")
        print(f"Максимальная скорость: {sol.y[2][-1]:.2f} м/с")
        print(f"Максимальная высота: {sol.y[3][-1]/1000:.2f} км")
        print(f"Угол траектории в конце: {sol.y[1][-1]*57.3:.2f}°")
        
        if attack_list:
            print(f"Максимальная атака: {max(attack_list):.2f}°")
            print(f"Атака в конце: {attack_list[-1]:.2f}°")
        
        # Сохранение результатов
        output()
        
        print(f"\nРасчет завершен успешно!")
        print(f"Сохрано {len(time_list)} точек данных")
        
    except Exception as e:
        print(f"\n❌ ОШИБКА ВО ВРЕМЯ РАСЧЕТА: {e}")
        import traceback
        traceback.print_exc()

# Запуск основной функции
if __name__ == "__main__":
    main()