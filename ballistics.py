import attack as att
import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp
import matplotlib.pyplot as plt
import math
import constants

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")

k1a = 4.1
k2a = 0.15

attack_list = []
time_list = []

class ballistics:
    def __init__(self, N, Y, vel, alt):
        self.N = N
        self.Y = Y
        self.vel = vel
        self.alt = alt
        
        time_sep = parser.get_work_time()
        self.alpha = att.alpha(k1a, k2a, time_sep[0], False)

        self.G = aero.UnionStream()
        self.G.set_elnumber(parser.get_block_number())
        self.G.set_diameter(parser.get_diameters())
        self.G.set_length(parser.get_part_length())

        self.thrust = 0
        self.mass = 0
        self.attack = 0

        
        self.atm = atmo.atmosphere(self.alt)

    def update_params(self, time):
        self.thrust =  parser.get_thrust_from_time(time)
        self.mass   =  parser.get_mass_from_time(time)
        self.attack = self.alpha.calculate_alpha(self.vel, time)*math.pi/180
        self.G.calculate_CXY(self.vel, self.alt, self.attack)

        self.atm = atmo.atmosphere(self.alt)

        attack_list.append(self.attack*180/math.pi)
        time_list.append(time)

        time += parser.interstep

    def delta_velocity(self, time):
        self.update_params(time)
        F_P = self.thrust * math.cos(self.attack)
        F_X = self.G.CX*self.atm.get_density()*parser.maximum_area*self.vel**2/2
        return F_P/self.mass-F_X/self.mass - self.atm.get_AOG() * math.sin(self.Y)

    def delta_trajangle(self, time):
        self.update_params(time)
        F_P = self.thrust * math.sin(self.attack)
        F_Y = self.G.CY*self.attack*self.atm.get_density()*parser.maximum_area*self.vel**2/2
        F_G = ((self.atm.get_AOG() * math.cos(self.Y)) * (1-self.vel**2/2)/(self.atm.get_AOG()*(constants.earth_radius + self.alt)))
        return (F_P + F_Y)/(self.mass*self.vel) - F_G/self.vel
    
    def delta_polar(self):
        return (self.vel/(constants.earth_radius + self.alt))*math.cos(self.Y)
    
    def delta_altitude(self):
        return self.vel * math.sin(self.Y)



from scipy.integrate import solve_ivp

v_max = 0


def system(t, vars):
    n, y, v, h = vars
    b = ballistics(n, y, v, h)
    dn = b.delta_polar()
    dy = b.delta_trajangle(t)
    dv = b.delta_velocity(t)
    dh = b.delta_altitude()
    return [dn, dy, dv, dh]

def event_stop_velocity(t, vars):
    v = vars[2]
    global v_max
    if (v > v_max): 
        v_max = v
    return v - 0.9

event_stop_velocity.terminal = True
event_stop_velocity.direction = -1

if __name__ == "__main__":
    ft = parser.get_full_time()
    h = parser.interstep
    t_span = (0, 120)
    y0 = [1, math.pi/2, 1, 1]
    
    sol = solve_ivp(system, t_span, y0, method='RK45', max_step=h, events=event_stop_velocity)
    print("Max velocity: ", v_max)

    plt.figure(figsize=(10,12))

    plt.subplot(4,1,1)
    plt.plot(sol.t, sol.y[2], label='Скорость v(t)', color='blue')
    plt.xlabel('Время, с')
    plt.ylabel('Скорость, м/с')
    plt.title('Скорость по времени')
    plt.legend()
    plt.grid(True)

    plt.subplot(4,1,2)
    plt.plot(time_list, attack_list, label='Угол атаки α(t)', color='red')
    plt.xlabel('Время, с')
    plt.ylabel('Угол атаки, градусы')
    plt.title('Угол атаки по времени')
    plt.legend()
    plt.grid(True)

    plt.subplot(4,1,3)
    plt.plot(sol.t, sol.y[3]/1000, label='Высота alt(t)', color='green')
    plt.xlabel('Время, с')
    plt.ylabel('Высота, км')
    plt.title('Высота по времени')
    plt.legend()
    plt.grid(True)

    plt.subplot(4,1,4)
    plt.plot(sol.t, sol.y[1]*57.3, label='Угол, град', color='blue')
    plt.xlabel('Время, с')
    plt.ylabel('Угол, град')
    plt.title('Угол по времени')
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()
