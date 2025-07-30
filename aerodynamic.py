import math

# Простая структура элемента ракеты
class Element:
    def __init__(self, upper_diameter=0, lower_diameter=0, length=0):
        self.upper_diameter = upper_diameter
        self.lower_diameter = lower_diameter
        self.length = length
        self.ratio = 0
        self.round_area = 0
        self.upper_area = 0
        self.lower_area = 0
        self.C_fric = 0
        self.C_pres = 0
        self.C_ind = 0
        self.CX = 0
        self.CY = 0

class Geometry:
    PI = math.pi

    def __init__(self):
        self.elements = []
        self.full_length = 0
        self.full_round_area = 0
        self.full_ratio = 0
        self.midel_area = 0

    def set_elements(self, diameters, lengths):
        n = len(diameters)
        self.elements = []
        for i in range(n):
            upper_d = diameters[i-1] if i > 0 else 0
            lower_d = diameters[i]
            length = lengths[i]
            elem = Element(upper_d, lower_d, length)
            self.elements.append(elem)
        self.pre_calculations()

    def pre_calculations(self):
        self.full_length = 0
        self.full_round_area = 0
        for i, e in enumerate(self.elements):
            if i == 0:
                e.base_line = 2 * self.PI * math.sqrt(e.length**2 + (e.upper_diameter/2)**2)
                e.round_area = 2 * self.PI * e.length * e.upper_diameter / 2 + self.PI * e.length**2
            else:
                e.base_line = math.sqrt(e.length**2 + ((e.lower_diameter - e.upper_diameter)**2)/4)
                e.round_area = self.PI * (e.upper_diameter + e.lower_diameter) * e.base_line / 2

            e.upper_area = self.PI * (e.upper_diameter/2)**2
            e.lower_area = self.PI * (e.lower_diameter/2)**2

            if e.upper_diameter < e.lower_diameter and abs(e.upper_diameter) > 0.1:
                e.ratio = (e.length + e.length / (e.lower_diameter / e.upper_diameter - 1)) / e.lower_diameter
            else:
                e.ratio = e.length / e.lower_diameter

            self.full_length += e.length
            self.full_round_area += e.round_area

        self.full_ratio = self.full_length / self.elements[-1].upper_diameter
        self.midel_area = self.elements[-1].lower_area

class Aerodynamics(Geometry):
    def __init__(self):
        super().__init__()

    def sqr(self, x):
        return x * x

    def rad(self, x):
        return x / 57.3

    def friction_calc(self, Mach, speed, nu):
        Re = speed * Mach * self.full_length / nu
        if Re <= 485000:
            cif = 2.656 / math.sqrt(Re)
            num = pow(1 + 0.1 * pow(Mach, 0.1), -0.125)
        elif Re < 1e7:
            n = 5 + (1.3 + 0.6 * Mach * (1 - 0.25 * pow(Mach, 2))) * \
                math.sqrt(1 - pow(math.log10((8e-6 / self.full_length * Re) - 1) / (2.2 + 0.08 * pow(Mach, 2) / (1 + 0.312 * pow(Mach, 2))), 2))
            x_t = min(pow(10, n) / Re, self.elements[0].length / self.full_length)
            if x_t >= 1:
                cif = 0.91 / pow(math.log10(Re), 2.58) * pow(1 - x_t + 40 * pow(x_t, 0.625) / pow(Re, 0.375), 0.8)
            else:
                cif = 2.656 / math.sqrt(Re)
            num = pow(1 + 0.1 * pow(Mach, 0.1), -2/3)
        else:
            cif = 0.91 / pow(math.log10(Re), 2.58)
            num = pow(1 + 0.1 * pow(Mach, 0.1), -2/3)

        area_ratio = self.full_round_area / self.midel_area
        for e in self.elements:
            e.C_fric = e.round_area * cif * num / 2

        C_fric = area_ratio * cif * num / 2
        return C_fric

    def pressure_calc(self, Mach):
        # Упрощённый пример: фиксированные коэффициенты для головы и конуса
        C_head = 0.035  # примерное значение
        C_cone = 0.02   # примерное значение

        C_pres = C_head + C_cone
        for e in self.elements[1:]:
            e.C_pres = C_cone * (1 - e.upper_area / e.lower_area)
            C_pres += e.C_pres
        return C_pres

    def lift_calc(self, Mach):
        # Упрощённый пример подъёмной силы
        C_lift = 0
        for e in self.elements:
            C_lift += 0.01 * Mach * e.ratio  # примерная зависимость
            e.CY = 0.01 * Mach * e.ratio
        return C_lift

    def calculate_CX(self, angle_deg, Mach, speed, nu):
        angle_rad = self.rad(angle_deg)
        Cx_fric = self.friction_calc(Mach, speed, nu)
        Cx_pres = self.pressure_calc(Mach)
        CY = self.lift_calc(Mach)
        E = 0.01  # упрощённое индуцированное сопротивление

        CX = Cx_fric + Cx_pres + CY + E
        CY = (CY + E) * angle_rad - self.rad(CY + E)
        CY *= angle_rad
        return CX, CY

# Пример использования
if __name__ == "__main__":
    # Задаём геометрию ракеты
    diameters = [3.7, 3.7, 4.1, 4.1]  # в метрах
    lengths = [7, 7, 2, 14]           # в метрах

    aero = Aerodynamics()
    aero.set_elements(diameters, lengths)

    # Параметры
    mach_numbers = [0.1, 0.3, 0.5, 0.7, 0.9]
    attack_angle_deg = 5
    speed = 340  # скорость звука, м/с (пример)
    nu = 1.5e-5  # кинематическая вязкость, м²/с (пример)

    for Mach in mach_numbers:
        CX, CY = aero.calculate_CX(attack_angle_deg, Mach, speed, nu)
        print(f"Mach={Mach:.2f}, CX={CX:.4f}, CY={CY:.4f}")
