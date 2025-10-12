import json
import attack
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

def calculate_static(mass_, shoulder):
    return 0.5 * mass_ * shoulder

def calculate_inertia(mass_, shoulder, shoulder_diff, diameter):
    return 0.25 * mass_ * (math.pow(shoulder, 2) + 0.333 * math.pow(shoulder_diff, 2) + math.pow((diameter/2), 2))

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
        self.diameters = r_data["diameters"]
        self.thrust = r_data["thrust"]
        self.attack_coefs = r_data["attack_coefs"]

        self.propellant_mass = []
        self.delta_mass = []
        self.stage_mass = []
        self.structural_mass = []
        self.delta_level_ox=[]
        self.delta_level_fu=[]
        self.delta_mass_ox=[]
        self.delta_mass_fu=[]
        self.work_time=[]
        self.full_time=[]

        self.thrust_vector = []
        self.mass_vector = []
        self.propellant_mass_vector = []
        self.time_vector = []
        self.mass_ox=[]
        self.mass_fu=[]

        # Added lists for static components
        self.static_ox = []
        self.static_str = []
        self.static_fu = []

        self.inertia_ox = []
        self.inertia_str = []
        self.inertia_fu = []

        mass_ox_t = []
        mass_fu_t = []

        time=0
        mass_t = self.full_mass
        thrust_t = self.thrust[0]

        for k in range(self.block_number):
            self.propellant_mass.append(self.block_mass[k]*self.structural_values[k]/(self.structural_values[k] + 1))
            self.mass_ox.append(self.propellant_mass[k]*self.components_ratio/(self.components_ratio+1))
            self.mass_fu.append(self.propellant_mass[k]*1/(self.components_ratio+1))
            self.delta_mass.append(self.thrust[k]/self.exhaust_velocity[k])
            self.work_time.append(self.propellant_mass[k]/self.delta_mass[k])

            self.stage_mass.append(self.full_mass - sum(self.block_mass[:k]))

            self.structural_mass.append(self.block_mass[k]-self.propellant_mass[k])

            self.delta_mass_ox.append(self.delta_mass[k]*self.components_ratio/(self.components_ratio+1))
            self.delta_mass_fu.append(self.delta_mass[k]*1/(self.components_ratio+1))

            self.delta_level_ox.append(self.delta_mass_ox[k]/self.oxidizer_density/self.maximum_area)
            self.delta_level_fu.append(self.delta_mass_fu[k]/self.fuel_density/self.maximum_area)

            mass_ox_t.append(self.mass_ox[k])
            mass_fu_t.append(self.mass_fu[k])

            shoulder_ox_t=(self.sector_index_ox[k][0] + self.sector_index_ox[k][-1])
            shoulder_fu_t=(self.sector_index_fu[k][0] + self.sector_index_fu[k][-1])
            shoulder_str =(self.sector_index_ox[k][0] + self.sector_index_fu[k][-1])

            shoulder_ox_diff_t  = abs(self.sector_index_ox[k][0] - self.sector_index_ox[k][-1])
            shoulder_fu_diff_t  = abs(self.sector_index_fu[k][0] - self.sector_index_fu[k][-1])
            shoulder_str_diff_t = abs(self.sector_index_ox[k][0] - self.sector_index_fu[k][-1])

            static_ox_t = calculate_static(mass_ox_t[k], shoulder_ox_t)
            static_fu_t = calculate_static(mass_fu_t[k], shoulder_fu_t)
            static_str  = calculate_static(self.structural_mass[k], shoulder_str)

            inertia_ox_t = calculate_inertia(mass_ox_t[k], shoulder_ox_t, shoulder_ox_diff_t, self.max_diameter)
            inertia_fu_t = calculate_inertia(mass_fu_t[k], shoulder_fu_t, shoulder_fu_diff_t, self.max_diameter)
            inertia_str  = calculate_inertia(self.structural_mass[k], shoulder_str, shoulder_str_diff_t, self.max_diameter)

            # Append to lists instead of accumulating a single static variable
            self.static_ox.append(static_ox_t)
            self.static_str.append(static_str)
            self.static_fu.append(static_fu_t)

            self.inertia_ox.append(inertia_ox_t)
            self.inertia_str.append(inertia_str)
            self.inertia_fu.append(inertia_fu_t)
            
            while mass_t >(self.stage_mass[k]-self.propellant_mass[k]):
                self.mass_vector.append(mass_t )
                self.time_vector.append(time)
                self.thrust_vector.append(thrust_t)
                mass_t -=self.delta_mass[k]*self.interstep
                time+=self.interstep
            mass_t -=self.structural_mass[k]
            if k+1>=(self.block_number):
                thrust_t = 0
            else:
                thrust_t = self.thrust[k+1]
        else:
            self.mass_vector.append(mass_t)
            self.time_vector.append(time)
            self.thrust_vector.append(thrust_t)
        
        # Calculate static as cumulative sums from k to end for each k
        self.full_time = sum(self.work_time)  
    
        self.static = []
        self.inertia = []
        for k in range(self.block_number):
            self.static.append(sum(self.static_ox[k:]) + sum(self.static_str[k:]) + sum(self.static_fu[k:]))
            print(self.static[k])
            self.inertia.append(sum(self.inertia_ox[k:]) + sum(self.static_str[k:]) + sum(self.static_fu[k:]))
            print(self.inertia[k])

        for m in range(len(self.time_vector)):
            for k in range(self.block_number):
                self.static[k] = 0

        
        self.alpha = attack.alpha(self.attack_coefs[0],
                                  self.attack_coefs[1],
                                  self.work_time[0],
                                  False)
    # rocket parameters
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
    def get_full_time(self):
        return self.full_time
    def get_diameters(self):
        return self.diameters
    def get_part_length(self):
        res = []
        res.append(self.head_length)
        for k in range(len(self.block_length)):
            res.append(self.block_length[k])
        return res
    def get_thrust(self):
        return self.thrust
    
    # flight parameters
    def vector_time(self):
        return self.time_vector
    def vector_mass(self):
        return self.mass_vector
    
    def vector_thrust(self):
        return self.thrust_vector
            
    def get_mass_from_time(self, time):
        for k in range(len(self.time_vector)):
            if abs(self.time_vector[k]-time)<self.interstep:
                return self.mass_vector[k]
    
    def get_thrust_from_time(self, time):
        for k in range(len(self.time_vector)):
            if abs(self.time_vector[k]-time)<self.interstep:
                return self.thrust_vector[k]

    def get_propellant_from_time(self, time):
        for k in range(len(self.time_vector)):
            if abs(self.time_vector[k]-time)<self.interstep:
                return self.thrust_vector[k]

    def get_attack(self, vel, time):
        return self.alpha.calculate_alpha(vel, time)