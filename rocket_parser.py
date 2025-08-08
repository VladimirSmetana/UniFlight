import json
import path
import math
import constants

def read_propellant_density(propellant_type):
    if propellant_type == "kerosene":
        return constants.density.kerosene.value
    elif propellant_type == "liquid_oxygen":
        return constants.density.liquid_oxygen.value
    elif propellant_type == "tetroxide":
        return constants.density.tetroxide.value
    elif propellant_type == "heptyl":
        return constants.density.heptyl.value
    else: return 0

class rocket_parser:
    def __init__(self, filename):
        with open(filename, 'r') as r_file:
            r_data = json.load(r_file)
        
        self.max_diameter = r_data["maximum_diameter"]
        self.head_length = r_data["head_length"]
        self.block_length = r_data["block_length"]
        self.element_length = r_data["rocket_sections"]
        self.rocket_length = self.head_length + sum(self.block_length)
        if abs(self.rocket_length - sum(self.element_length)) > 0.001:
            print("Rocket length Error: " + str(self.head_length + sum(self.block_length)) + "!=" + str(sum(self.element_length)))

        self.maximum_area = (math.pi*(self.max_diameter**2))/4
        self.structural_values = r_data["structural_values"]
        self.payload = r_data["payload_mass"]
        self.block_mass = r_data["block_mass"]
        self.full_mass = self.payload + sum(self.block_mass)

        self.block_number = r_data["block_number"]
        self.thrust = r_data["thrust"]
        self.exhaust_velocity = r_data["exhaust_velocity"]
        self.components_ratio = r_data["components_ratio"]
        self.sector_index = r_data["sector_index"]
        self.sector_index_ox = self.sector_index["ox"]
        self.sector_index_fu = self.sector_index["fu"]
        self.fuel_density = read_propellant_density(r_data["fuel"])
        self.oxidizer_density = read_propellant_density(r_data["oxidizer"])
        self.interstep = r_data["integration_step"]

        self.propellant_mass = []
        self.delta_mass = []
        self.stage_mass = []
        self.structural_mass = []
        self.delta_level_ox=[]
        self.delta_level_fu=[]
        self.delta_mass_ox=[]
        self.delta_mass_fu=[]
        self.work_time=[]

        self.stage_mass.append(self.full_mass)

        for k in range(self.block_number):
            self.propellant_mass.append(self.block_mass[k]*self.structural_values[k]/(self.structural_values[k] + 1))
            self.delta_mass.append(self.thrust[k]/self.exhaust_velocity[k])
            self.work_time.append(self.propellant_mass[k]/self.delta_mass[k])

            self.stage_mass.append(self.full_mass-self.block_mass[k])

            self.structural_mass.append(self.block_mass[k]-self.propellant_mass[k])

            self.delta_mass_ox.append(self.delta_mass[k]*self.components_ratio/(self.components_ratio+1))
            self.delta_mass_fu.append(self.delta_mass[k]*1/(self.components_ratio+1))

            self.delta_level_ox.append(self.delta_mass_ox[k]/self.oxidizer_density/self.maximum_area)
            self.delta_level_fu.append(self.delta_mass_fu[k]/self.fuel_density/self.maximum_area)

    def get_block_number(self):
        return self.block_number
    def get_rocket_length(self):
        return self.rocket_length
    def get_structural_values(self):
        return self.structural_values
    def get_structural_mass(self):
        return self.structural_mass
    def get_payload(self):
        return self.payload
    def get_full_mass(self):
        return self.full_mass
    def get_propellant_mass(self):
        return self.propellant_mass
    def get_delta_mass(self):
        return self.delta_mass
    def get_stage_mass(self):
        return self.stage_mass
    def get_delta_level_ox(self):
        return self.delta_level_ox    
    def get_delta_level_fu(self):
        return self.delta_level_fu 
    def get_delta_mass_ox(self):
        return self.delta_mass_ox    
    def get_delta_mass_fu(self):
        return self.delta_mass_fu 
    def get_interstep(self):
        return self.interstep
    def get_sector_index_ox(self):
        return self.sector_index_ox
    def get_sector_index_fu(self):
        return self.sector_index_fu
    def get_work_time(self):
        return self.work_time