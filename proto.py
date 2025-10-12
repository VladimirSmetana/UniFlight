import atmosphere as atmo
import aerodynamics as aero
import path
import rocket_parser as rp
import matplotlib.pyplot as plt

parser = rp.rocket_parser(path.rocket_lib + "master_rocket.json")


full_mass = parser.get_full_mass()
delta_mass = parser.get_delta_mass()
stage_mass = parser.get_stage_mass()
propellant_mass = parser.get_propellant_mass()
structural_mass = parser.get_structural_mass()
block_number = parser.get_block_number()
time_points = parser.get_work_time()


plt.plot(parser.vector_time(), parser.static_time)
plt.title('Расчет массы РН по времени полета', fontsize=16)
plt.xlabel('Время полета, с', fontsize=14)
plt.ylabel('Масса, кг', fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()