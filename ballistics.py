import attack as att
import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp
import matplotlib.pyplot as plt
import math
import constants

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")


class ballistics:
    def __init__(self, N, Y, vel):
        self.N = N
        self.Y = Y
        self.vel = vel

        self.alt = self.vel * parser.interstep * math.sin(self.Y)

        self.alpha = att.alpha(1, 1, 2, 0)

        self.G = aero.UnionStream()
        self.G.set_elnumber(parser.get_block_number())
        self.G.set_diameter(parser.get_diameters())
        self.G.set_length(parser.get_part_length())

        self.thrust = 0
        self.mass = 0
        self.attack = 0
        self.atm = atmo.atmosphere(0)

    def update_params(self, time):
        print("getting mass ps")
        self.thrust =  parser.get_thrust_from_time(time)
        self.mass   =  parser.get_mass_from_time(time)
        print("getting attack angle")
        self.attack = self.alpha.calculate_alpha(self.vel, time)
        print("getting cx")
        self.G.calculate_CXY(self.vel, self.alt, self.attack)
        print("getting atmosphere ps")
        self.atm = atmo.atmosphere(self.alt)

    def delta_velocity(self, time):
        self.update_params(time)
        F_P = self.thrust * math.cos(self.attack)
        F_X = self.G.CX*self.atm.get_density()*parser.maximum_area*self.vel**2/2
        return F_P/self.mass-F_X/self.mass - self.atm.get_AOG() * math.sin(self.Y)

    def delta_trajangle(self, time):
        self.update_params(time)
        F_P = self.thrust * math.sin(self.attack)
        F_Y = self.G.CY*self.atm.get_density()*parser.maximum_area*self.vel**2/2
        F_G = ((self.atm.get_AOG() * math.cos(self.Y)) * (1-self.vel**2/2)/(self.atm.get_AOG()*(constants.earth_radius + self.alt)))
        return (F_P + F_Y)/(self.mass*self.vel) - F_G/self.vel
    
    def delta_polar(self):
        return (self.vel/(constants.earth_radius + self.alt))*math.cos(self.Y)


if __name__ == "__main__":

    n = 0
    y = 0
    v = 0
    t = 0

    h = parser.interstep

    while t < 10:
        b = ballistics(n, y, v)
        v = h * b.delta_velocity(t)
        y = h * b.delta_trajangle(t)
        n = h * b.delta_polar()
        t += h