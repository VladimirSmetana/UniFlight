import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp
import matplotlib.pyplot as plt

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")


full_mass = parser.get_full_mass()
delta_mass = parser.get_delta_mass()
stage_mass = parser.get_stage_mass()
fuel_mass = parser.get_fuel_mass()
structural_mass = parser.get_structural_mass()
block_number = parser.get_block_number()

h=0.01
time=0
mass_vector = []
time_vector = []

for k in range(block_number):
    while full_mass>(stage_mass[k]-fuel_mass[k]):
        mass_vector.append(full_mass)
        time_vector.append(time)
        full_mass-=delta_mass[k]*h
        time+=h
    full_mass-=structural_mass[k]
else:
    mass_vector.append(full_mass)
    time_vector.append(time)

print(full_mass)

plt.plot(time_vector, mass_vector)
plt.show()