import math

class alpha:
    def __init__(self, k1: float, k2: float, septime: float, IsRichOrbit: bool):
        self.k1 = k1
        self.k2 = k2
        self.septime = septime
        self.IsRichOrbit = IsRichOrbit

    def calculate_alpha(self, velocity: float, time: float) -> float:
        ans = 0.0
        z = 0.0
        che = 0.0

        if 50 < velocity < 270 and time <= self.septime:
            z = math.pi * (velocity - 50)
            che = (velocity - 50) + 0.25 * (270 - velocity)
            ans = - self.k1 * (math.sin(z / che) ** 2)
        elif time >= self.septime:
            if not self.IsRichOrbit:
                ans = self.k2 * (time - self.septime)**(1/2)
            else:
                if time - self.septime < 60:
                    z = math.pi * (time - self.septime)
                    che = (time - self.septime) + 0.25 * (self.septime + 60 - time)
                    ans = 90 * (math.sin(z / che) ** 2)
        return ans

# Пример использования
if __name__ == "__main__":
    a = alpha(1.0, 0.5, 100.0, True)
    print(a.calculate_alpha(120.0, 95.0))