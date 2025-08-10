import enum

acceleration_of_gravity = 9.81
earth_radius = 6371000
lamb = (4.73, 7.853, 10.996, 14.137, 17.279)

class density(enum.Enum):
    liquid_oxygen = 1100.0
    kerosene  = 440.0
    tetroxide = 1450.0
    heptyl = 790.0