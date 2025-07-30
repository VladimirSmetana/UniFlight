import math
import matplotlib.pyplot as plt
import atmosphere

class Element:
    def __init__(self):
        self.upper_diameter = 0.0
        self.lower_diameter = 0.0
        self.elem_length = 0.0
        self.virtual_length = 0.0
        self.ratio = 0.0
        self.round_area = 0.0
        self.upper_area = 0.0
        self.lower_area = 0.0
        self.base_line = 0.0

        self.C_fric = 0.0
        self.C_pres = 0.0
        self.C_ind = 0.0
        self.CX = 0.0
        self.CY = 0.0

class Geometry:
    PI = math.pi

    def __init__(self):
        self.elem = []
        self.full_length = 0.0
        self.full_round_area = 0.0
        self.full_ratio = 0.0
        self.midel_diameter = 0.0
        self.midel_area = 0.0
        self.cif = 0.0
        self.num = 0.0

    def sqr(self, x):
        return x * x

    def rad(self, x):
        return x / 57.3

    def set_elnumber(self, n):
        self.elem = [Element() for _ in range(n)]

    def set_length(self, lengths):
        for i, l in enumerate(lengths):
            self.elem[i].elem_length = l
        self.pre_calculations()

    def set_diameter(self, diameters):
        for i, d in enumerate(diameters):
            if i == 0:
                self.elem[i].upper_diameter = 0.0
                self.elem[i].lower_diameter = d
            else:
                self.elem[i].upper_diameter = diameters[i - 1]
                self.elem[i].lower_diameter = d

    def pre_calculations(self):
        self.full_length = 0.0
        self.full_round_area = 0.0
        for i, e in enumerate(self.elem):
            if i == 0:
                e.base_line = 2 * self.PI * math.sqrt(e.elem_length ** 2 + (e.upper_diameter / 2) ** 2)
                e.round_area = 2 * self.PI * e.elem_length * e.upper_diameter / 2 + self.PI * (e.elem_length ** 2)
            else:
                e.base_line = math.sqrt(e.elem_length ** 2 + ((e.lower_diameter - e.upper_diameter) ** 2) / 4)
                e.round_area = self.PI * (e.upper_diameter + e.lower_diameter) * e.base_line / 2

            e.upper_area = self.PI * (e.upper_diameter / 2) ** 2
            e.lower_area = self.PI * (e.lower_diameter / 2) ** 2

            try:
                e.virtual_length = e.elem_length + e.elem_length / (e.lower_diameter / e.upper_diameter - 1)
            except ZeroDivisionError:
                e.virtual_length = e.elem_length

            self.full_length += e.elem_length
            self.full_round_area += e.round_area

            if e.upper_diameter < e.lower_diameter and abs(e.upper_diameter) > 0.1:
                e.ratio = e.virtual_length / e.lower_diameter
            else:
                if e.lower_diameter != 0:
                    e.ratio = e.elem_length / e.lower_diameter
                else:
                    e.ratio = 0

        if self.elem[-1].upper_diameter != 0:
            self.full_ratio = self.full_length / self.elem[-1].upper_diameter
        else:
            self.full_ratio = 0

        self.midel_diameter = self.elem[-1].upper_diameter
        self.midel_area = self.elem[-1].lower_area

class Friction(Geometry):
    def __init__(self):
        super().__init__()
        self.h_s = 8e-6
        self.area_ratio = 0.0
        self.Re = 0.0
        self.n = 0.0
        self.x_t = 0.0
        self.C_fric = 0.0

    def stream_calc(self, Re, Mach):
        if Re <= 485000:
            self.cif = 2.656 / math.sqrt(Re)
            self.num = pow(1 + 0.1 * pow(Mach, 0.1), -0.125)
        elif Re < 10000000:
            try:
                log_arg = (self.h_s / self.full_length * Re) - 1
                denominator = 2.2 + 0.08 * pow(Mach, 2) / (1 + 0.312 * pow(Mach, 2))
                self.n = 5 + (1.3 + 0.6 * Mach * (1 - 0.25 * pow(Mach, 2))) * math.sqrt(
                    1 - pow(math.log10(log_arg) / denominator, 2))
            except (ValueError, ZeroDivisionError):
                self.n = 5

            self.x_t = min(pow(10, self.n) / Re, self.elem[0].elem_length / self.full_length)
            if self.x_t >= 1:
                self.cif = 0.91 / pow(math.log10(Re), 2.58) * pow(1 - self.x_t + 40 * pow(self.x_t, 0.625) / pow(Re, 0.375), 0.8)
            else:
                self.cif = 2.656 / math.sqrt(Re)
            self.num = pow(1 + 0.1 * pow(Mach, 0.1), -2 / 3)
        else:
            self.cif = 0.91 / pow(math.log10(Re), 2.58)
            self.num = pow(1 + 0.1 * pow(Mach, 0.1), -2 / 3)

    def fricalc(self, Mach, SS, nu):
        if self.midel_area == 0:
            self.area_ratio = 0
        else:
            self.area_ratio = self.full_round_area / self.midel_area
        self.Re = SS * Mach * self.full_length / nu if nu != 0 else 0
        self.stream_calc(self.Re, Mach)

        for e in self.elem:
            e.C_fric = e.round_area * self.cif * self.num / 2

        self.C_fric = self.area_ratio * self.cif * self.num / 2
        return self.C_fric

class Pressure(Geometry):
    def read_pressure_file(self, filename, rows, cols):
        data = []
        with open(filename, 'r') as f:
            for _ in range(rows):
                line = f.readline()
                if not line:
                    break
                parts = line.strip().split()
                if len(parts) < cols:
                    raise ValueError(f"File {filename} line {_+1} has fewer than {cols} columns")
                data.append([float(x) for x in parts[:cols]])
        transposed = list(map(list, zip(*data)))
        return transposed

    def interpolate_Mach(self, Mach, data):
        Mach_v = data[0]
        values = data[1]
        for i in range(1, len(Mach_v)):
            if Mach >= Mach_v[i - 1] and Mach < Mach_v[i]:
                return values[i - 1] + (Mach - Mach_v[i - 1]) * (values[i] - values[i - 1]) / (Mach_v[i] - Mach_v[i - 1])
        return values[-1]

    def select_ratio_data_pressure(self, ratio, data):
        # data: [Mach_v, H_0, H_025, H_05, H_1, H_2, H_25, H_3, H_4, H_5]
        if ratio < 0.25:
            return [data[0], data[1]]
        elif ratio < 0.5:
            return [data[0], data[2]]
        elif ratio < 1:
            return [data[0], data[3]]
        elif ratio < 2:
            return [data[0], data[4]]
        elif ratio < 2.5:
            return [data[0], data[5]]
        elif ratio < 3:
            return [data[0], data[7]]
        elif ratio < 4:
            return [data[0], data[8]]
        else:
            return [data[0], data[9]]

    def select_ratio_data_triangle(self, ratio, data):
        if ratio >= 1.5 and ratio < 2:
            return [data[0], data[1]]
        elif ratio >= 2 and ratio < 2.5:
            return [data[0], data[2]]
        elif ratio >= 2.5 and ratio < 3:
            return [data[0], data[3]]
        elif ratio >= 3 and ratio < 4:
            return [data[0], data[4]]
        elif ratio >= 4:
            return [data[0], data[5]]
        else:
            return [data[0], data[1]]

    def bottom_pres(self, Mach):
        if Mach < 1:
            if self.cif * self.num * self.full_ratio == 0:
                return 0
            return 0.0155 / math.sqrt(self.cif * self.num * self.full_ratio)
        else:
            data = self.read_pressure_file("HeadPressure.txt", 10, 10)
            # 7-й столбец (индекс 6) - как в исходнике
            return self.interpolate_Mach(Mach, [data[0], data[6]])

    def head_Cpres(self, Mach):
        data = self.read_pressure_file("HeadPressure.txt", 10, 10)
        ratio = self.elem[0].ratio
        H_current = self.select_ratio_data_pressure(ratio, data)
        return self.interpolate_Mach(Mach, H_current)

    def triangle_Cpres(self, Mach, ratio):
        data = self.read_pressure_file("TrianglePressure.txt", 10, 7)
        H_current = self.select_ratio_data_triangle(ratio, data)
        return self.interpolate_Mach(Mach, H_current)

    def prescalc(self, Mach):
        res = 0.0
        for i in range(1, len(self.elem)):
            self.elem[i].C_pres = self.triangle_Cpres(Mach, self.elem[i].ratio) * (1 - self.elem[i].upper_area / self.elem[i].lower_area)
            res += self.elem[i].C_pres
        return res + self.head_Cpres(Mach) + self.bottom_pres(Mach)

class Inductance(Geometry):
    def read_pressure_file(self, filename, rows, cols):
        data = []
        with open(filename, 'r') as f:
            for _ in range(rows):
                line = f.readline()
                if not line:
                    break
                parts = line.strip().split()
                if len(parts) < cols:
                    raise ValueError(f"File {filename} line {_+1} has fewer than {cols} columns")
                data.append([float(x) for x in parts[:cols]])
        transposed = list(map(list, zip(*data)))
        return transposed

    def sqr(self, x):
        return x * x

    def rad(self, x):
        return x / 57.3

    def E_pressure(self, angle, Mach):
        N = 10
        data = self.read_pressure_file("EPressure.txt", N, 3)
        Mach_v = data[0]
        H_head = data[1]
        H_cone = data[2]

        if Mach < 1:
            Mach_val = -math.sqrt(1 - self.sqr(Mach)) / self.elem[0].ratio if self.elem[0].ratio != 0 else 0
        else:
            Mach_val = math.sqrt(self.sqr(Mach) - 1) / self.elem[0].ratio if self.elem[0].ratio != 0 else 0

        E_head = 0
        E_cone = 0
        for i in range(1, N):
            if Mach_val >= Mach_v[i - 1] and Mach_val < Mach_v[i]:
                E_head = H_head[i - 1] + (Mach_val - Mach_v[i - 1]) * (H_head[i] - H_head[i - 1]) / (Mach_v[i] - Mach_v[i - 1])
                E_cone = H_cone[i - 1] + (Mach_val - Mach_v[i - 1]) * (H_cone[i] - H_cone[i - 1]) / (Mach_v[i] - Mach_v[i - 1])

        self.elem[0].C_ind = (self.elem[0].CY + self.rad(2 * E_head)) * self.sqr(angle)
        for j in range(1, len(self.elem)):
            if self.elem[j].upper_diameter < self.elem[j].lower_diameter:
                ratio = self.elem[-1].upper_area if self.elem[-1].upper_area != 0 else 1
                self.elem[j].C_ind = (self.elem[j].CY * self.elem[j].upper_area / ratio + self.rad(2 * E_cone * self.elem[j].upper_area / ratio)) * self.sqr(angle)

        E = sum(e.C_ind for e in self.elem)
        return E

class LiftForce(Inductance):
    def sqr(self, x):
        return x * x

    def rad(self, x):
        return x / 57.3

    def head_lift(self, Mach):
        N = 9
        data = self.read_pressure_file("HeadNormal.txt", N, 6)
        Mah_v = data[0]
        H_0 = data[1]
        H_05 = data[2]
        H_1 = data[3]
        H_2 = data[4]
        H_4 = data[5]

        L_cyl = 0.0
        for j in range(1, len(self.elem)):
            if self.elem[j].upper_diameter < self.elem[j].lower_diameter:
                if self.elem[j].lower_diameter != 0:
                    L_cyl /= self.elem[j].lower_diameter
                break
            else:
                L_cyl += self.elem[j].elem_length

        ratio = L_cyl / self.elem[0].ratio if self.elem[0].ratio != 0 else 0

        H_current = []
        for i in range(N):
            print(i)
            if 0 <= ratio < 0.5:
                H_current.append(H_0[i])
            elif 0.5 <= ratio < 1:
                H_current.append(H_05[i])
            elif 1 <= ratio < 2:
                H_current.append(H_1[i])
            elif 2 <= ratio < 4:
                H_current.append(H_2[i])
            elif ratio >= 4:
                H_current.append(H_4[i])
            else:
                H_current.append(H_0[i])

        if Mach < 1:
            Mach_val = -math.sqrt(1 - self.sqr(Mach)) / self.elem[0].ratio if self.elem[0].ratio != 0 else 0
        else:
            Mach_val = math.sqrt(self.sqr(Mach) - 1) / self.elem[0].ratio if self.elem[0].ratio != 0 else 0

        C_head = 0.035
        for i in range(1, N):
            if Mach_val >= Mah_v[i - 1] and Mach_val < Mah_v[i]:
                C_head = H_current[i - 1] + (Mach_val - Mah_v[i - 1]) * (H_current[i] - H_current[i - 1]) / (Mah_v[i] - Mah_v[i - 1])

        self.elem[0].CY = C_head
        return C_head

    def free_triangle_lift(self, index):
        arg = self.elem[index].lower_diameter / 2 / (self.elem[index].virtual_length - self.elem[index].elem_length) if (self.elem[index].virtual_length - self.elem[index].elem_length) != 0 else 0
        Q = math.atan(arg)
        return (2 / 57.3) * self.sqr(math.cos(Q))

    def triangle_lift(self, Mach, ratio, index):
        N = 10
        data = self.read_pressure_file("TriangleNormal.txt", N, 6)
        Mah_v = data[0]
        H_0 = data[1]
        H_1 = data[2]
        H_2 = data[3]
        H_3 = data[4]
        H_4 = data[5]

        L_cyl = 0.0
        for j in range(index + 1, len(self.elem)):
            if self.elem[j].upper_diameter < self.elem[j].lower_diameter:
                if self.elem[j].lower_diameter != 0:
                    L_cyl /= self.elem[j].lower_diameter
                break
            else:
                L_cyl += self.elem[j].elem_length

        ratio_new = L_cyl / ratio if ratio != 0 else 0

        H_current = []
        for i in range(N):
            if 0 < ratio_new < 1:
                H_current.append(H_0[i])
            elif 1 <= ratio_new < 2:
                H_current.append(H_1[i])
            elif 2 <= ratio_new < 3:
                H_current.append(H_2[i])
            elif 3 <= ratio_new < 4:
                H_current.append(H_3[i])
            elif ratio_new >= 4:
                H_current.append(H_4[i])
            else:
                H_current.append(H_0[i])

        if Mach < 1:
            Mach_val = -math.sqrt(1 - self.sqr(Mach)) / ratio if ratio != 0 else 0
        else:
            Mach_val = math.sqrt(self.sqr(Mach) - 1) / ratio if ratio != 0 else 0

        C_head = 0.0
        for i in range(1, N):
            if Mach_val >= Mah_v[i - 1] and Mach_val < Mah_v[i]:
                C_head = H_current[i - 1] + (Mach_val - Mah_v[i - 1]) * (H_current[i] - H_current[i - 1]) / (Mah_v[i] - Mah_v[i - 1])

        return C_head

    def calculate_CY(self, Mach):
        return self.head_lift(Mach) + self.un_triangle_lift(Mach)

    def un_triangle_lift(self, Mach):
        res = 0.0
        for i in range(1, len(self.elem)):
            if self.elem[i].upper_area < self.elem[i].lower_area:
                big_rat = self.elem[i].ratio
                S_rat = self.elem[i].upper_area / self.elem[i].lower_area if self.elem[i].lower_area != 0 else 0
                self.elem[i].CY = self.triangle_lift(Mach, big_rat, i) - self.free_triangle_lift(i) * S_rat
                res += self.elem[i].CY * self.elem[i].upper_area / self.elem[-1].upper_area if self.elem[-1].upper_area != 0 else 0
        return res

class DragForce(Friction, Pressure):
    def calculate_CX(self, Mach, SS, nu):
        return self.fricalc(Mach, SS, nu) + self.prescalc(Mach)

class UnionStream(DragForce, LiftForce):
    def __init__(self):
        super().__init__()
        self.E = 0.0
        self.CX = 0.0
        self.CY = 0.0

    def calculate_CXY(self, angle, Mach, SS, nu):
        self.CX = self.calculate_CX(Mach, SS, nu)
        self.CY = self.calculate_CY(Mach)
        self.E = self.E_pressure(angle, Mach)
        self.CX += (self.CY + self.E)
        self.CY -= self.rad(self.CY + self.E)
        self.CY *= angle

def main():
    G = UnionStream()
    G.set_elnumber(4)
    G.set_diameter([3.7, 3.7, 4.1, 4.1])
    G.set_length([7, 7, 2, 14])

    arrayMach = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.3, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
    attack_angles = [2, 4, 6, 8, 10]  # градусы
    altitudes = [0, 10, 20, 40, 60]  # км

    plt.figure(figsize=(14, 6))

    for j in range(len(attack_angles)):
        A = atmosphere.Init(altitudes[j]*1000)
        angle_rad = G.rad(attack_angles[j])
        CX_list = []
        CY_list = []
        for mach in arrayMach:
            G.calculate_CXY(angle_rad, mach, A.get_SV(), A.get_dyn())
            CX_list.append(G.CX)
            CY_list.append(G.CY)

        plt.subplot(1, 2, 1)
        plt.plot(arrayMach, CX_list, label=f'α={attack_angles[j]}°, H={altitudes[j]}km')
        plt.xlabel('Mach')
        plt.ylabel('CX')
        plt.grid(True)
        plt.legend(fontsize=8)

        plt.subplot(1, 2, 2)
        plt.plot(arrayMach, CY_list, label=f'α={attack_angles[j]}°, H={altitudes[j]}km')
        plt.xlabel('Mach')
        plt.ylabel('CY')
        plt.grid(True)
        plt.legend(fontsize=8)

    plt.suptitle('Аэродинамические коэффициенты CX и CY')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

if __name__ == "__main__":
    main()
