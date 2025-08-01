import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")

rocket_length = parser.get_rocket_length()
structural_values = parser.get_structural_values()

print(parser.get_fuel_mass())