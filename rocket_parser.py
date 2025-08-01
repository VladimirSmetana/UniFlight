import json
import path


class rocket_parser:
    def __init__(self, filename):
        with open(filename, 'r') as r_file:
            r_data = json.load(r_file)
        
        self.head_length = r_data["head_length"]
        self.block_length = r_data["block_length"]
        self.element_length = r_data["rocket_sections"]
        self.rocket_length = self.head_length + sum(self.block_length)
        if abs(self.rocket_length - sum(self.element_length)) > 0.001:
            print("Rocket length Error: " + str(self.head_length + sum(self.block_length)) + "!=" + str(sum(self.element_length)))

        self.structural_values = r_data["structural_values"]
        self.payload = r_data["payload_mass"]
        self.block_mass = r_data["block_mass"]
        self.full_mass = self.payload + sum(self.block_mass)

        self.block_number = r_data["block_number"]
        self.fuel_mass = []
        for k in range(self.block_number):
            self.fuel_mass.append(self.block_mass[k]*self.structural_values[k]/(self.structural_values[k] + 1))

    def get_rocket_length(self):
        return self.rocket_length
    def get_structural_values(self):
        return self.structural_values
    def get_payload(self):
        return self.payload
    def get_full_mass(self):
        return self.full_mass
    def get_fuel_mass(self):
        return self.fuel_mass
        
