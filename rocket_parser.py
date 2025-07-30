import json



def parse_rocket(filename):
    with open(filename, 'r') as r_file:
        r_data = json.load(r_file)
        
        head_length = r_data["head_length"]
        block_length = r_data["block_length"]
        element_length = r_data["rocket_sections"]
        if abs(head_length + sum(block_length) - sum(element_length)) > 0.001:
            print("Rocket length Error: " + str(head_length + sum(block_length)) + "!=" + str(sum(element_length)))
        



if __name__ == "__main__":
    parse_rocket('master_rocket.json')