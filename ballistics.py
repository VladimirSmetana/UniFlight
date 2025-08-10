import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp
import matplotlib.pyplot as plt

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")


class classs(N, Y, vel)

def delta_velocity(time):
    F_P = P(t) * cos(attack(time, vel))
    F_X = CX(vel, alt, attack)*po(alt)*vel**2/2
    return F_P/mass(time)-F_X/mass(time) - g(alt) * sin(Y...)


    thrust = parser.get_thrust_from_time(time)
    return thrust
    # double equations::fdV(double vv, double ii)
    # {
    #     F_P = (PENG * cos((M_PI * alpha) / 180));
    #     F_X = CX * S * po * pow(vv, 2) / 2;
    #     return  F_P/m -  F_X/m - g * sin(ii);
    # }

 
    # double equations::fdY(double hh, double vv, double ii)
    # {
    #     F_P = (PENG * sin((M_PI * alpha) / 180));
    #     F_Y = (CY * ((M_PI * alpha) / 180) * S * (po * pow(vv, 2)) / 2);
    #     return (F_P+ F_Y)/ (m*vv)  - ((g  * cos(ii))) * (1-pow(vv,2)/(g*(constants::earth_radius+hh))) /vv;
    # }


    # double equations::fdN(double hh, double vv, double ii)
    # {
    #     return (vv /(constants::earth_radius+hh)) * cos(ii);
    # }

def update_velocity():
def update_anomaly_angle():
def update_polar_angle():

if __name__ == "__main__":
    print(delta_velocity(200))