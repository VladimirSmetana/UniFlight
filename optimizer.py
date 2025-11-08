import math
from scipy.optimize import minimize
import numpy as np
import attack as alpha

class AlphaOptimizer:
    def __init__(self, target_altitude=210000, target_velocity=7800, 
                 max_attack_angle=5, septime_range=(80, 150)):
        self.target_altitude = target_altitude
        self.target_velocity = target_velocity  # первая космическая ~7800 м/с
        self.max_attack_angle = max_attack_angle
        self.septime_range = septime_range
        
    def evaluate_trajectory(self, k1, k2, septime, IsRichOrbit):
        """Запускает баллистическое моделирование и возвращает оценку качества"""
        try:
            # Создаем alpha с текущими параметрами
            alpha_obj = alpha(k1, k2, septime, IsRichOrbit)
            
            # Запускаем баллистическое моделирование (упрощенная версия)
            final_conditions = self.simulate_trajectory(alpha_obj)
            
            # Вычисляем функцию стоимости
            cost = self.calculate_cost(final_conditions)
            return cost
            
        except Exception as e:
            return float('inf')
    
    def simulate_trajectory(self, alpha_obj):
        """Упрощенное моделирование траектории"""
        # Здесь должна быть интеграция с ballistics классом
        # Для примера - упрощенная версия
        
        # Имитация конечных условий полета
        # В реальности это будет вызов solve_ivp как в основном коде
        final_altitude = 0
        final_velocity = 0
        final_attack = 0
        
        # Заглушка - в реальности здесь полное моделирование
        return {
            'altitude': final_altitude,
            'velocity': final_velocity, 
            'attack_angle': final_attack
        }
    
    def calculate_cost(self, conditions):
        """Вычисляет функцию стоимости для оптимизации"""
        alt_error = abs(conditions['altitude'] - self.target_altitude) / 1000  # км
        vel_error = max(0, self.target_velocity - conditions['velocity']) / 1000  # км/с
        attack_error = abs(conditions['attack_angle'])  # градусы
        
        # Штрафуем за невыполнение условий
        cost = (alt_error * 10 +  # высота важнее
                vel_error * 100 +  # скорость очень важна
                attack_error * 5 +  # угол атаки важен
                (1 if conditions['velocity'] < self.target_velocity else 0) * 1000)  # большой штраф за недостаток скорости
        
        return cost
    
    def optimize_alpha(self, initial_guess=None):
        """Оптимизирует параметры alpha"""
        if initial_guess is None:
            initial_guess = [1.0, 0.5, 100.0, False]  # [k1, k2, septime, IsRichOrbit]
        
        # Границы параметров
        bounds = [
            (0.1, 10.0),    # k1
            (0.01, 5.0),    # k2  
            self.septime_range,  # septime
            (0, 1)          # IsRichOrbit как 0/1
        ]
        
        def objective(params):
            k1, k2, septime, is_rich_int = params
            is_rich = bool(is_rich_int)
            return self.evaluate_trajectory(k1, k2, septime, is_rich)
        
        # Запуск оптимизации
        result = minimize(objective, initial_guess, bounds=bounds, 
                         method='L-BFGS-B', options={'maxiter': 100})
        
        if result.success:
            optimized_params = result.x
            return {
                'k1': optimized_params[0],
                'k2': optimized_params[1], 
                'septime': optimized_params[2],
                'IsRichOrbit': bool(optimized_params[3]),
                'cost': result.fun
            }
        else:
            raise Exception("Оптимизация не удалась")

class AdaptiveAlpha(alpha):
    """Расширенный класс alpha с автоматической оптимизацией"""
    
    def __init__(self, k1: float = None, k2: float = None, septime: float = None, 
                 IsRichOrbit: bool = None, auto_optimize=True):
        
        if auto_optimize and (k1 is None or k2 is None or septime is None):
            # Автоматическая оптимизация параметров
            optimizer = AlphaOptimizer()
            optimized = optimizer.optimize_alpha()
            
            k1 = optimized['k1']
            k2 = optimized['k2'] 
            septime = optimized['septime']
            IsRichOrbit = optimized['IsRichOrbit']
            
            print(f"Автоматически оптимизированные параметры:")
            print(f"k1={k1:.3f}, k2={k2:.3f}, septime={septime:.1f}, IsRichOrbit={IsRichOrbit}")
        
        super().__init__(k1, k2, septime, IsRichOrbit)
        self.is_optimized = auto_optimize

# Модифицированный класс ballistics с проверкой условий
class AdaptiveBallistics(ballistics):
    def __init__(self, N, Y, vel, alt):
        super().__init__(N, Y, vel, alt)
        self.optimizer = AlphaOptimizer()
        
    def check_target_conditions(self, time):
        """Проверяет достижение целевых условий"""
        self.update_params(time)
        
        current_altitude = self.alt
        current_velocity = self.vel
        current_attack = abs(self.attack * 180 / math.pi)  # в градусах
        
        # Проверка условий
        altitude_ok = 200000 <= current_altitude <= 220000
        velocity_ok = current_velocity >= 7800  # первая космическая
        attack_ok = current_attack <= 5
        
        return altitude_ok and velocity_ok and attack_ok
    
    def optimize_if_needed(self, current_time):
        """Запускает оптимизацию если условия не выполняются"""
        if current_time > 100:  # Даем время на разгон
            if not self.check_target_conditions(current_time):
                print("Условия не выполняются, запуск оптимизации...")
                optimized_params = self.optimizer.optimize_alpha()
                
                # Обновляем параметры alpha
                self.attack = AdaptiveAlpha(
                    optimized_params['k1'],
                    optimized_params['k2'], 
                    optimized_params['septime'],
                    optimized_params['IsRichOrbit'],
                    auto_optimize=False
                )

# Использование в основном коде
if __name__ == "__main__":
    # Автоматическая оптимизация при создании
    optimized_alpha = AdaptiveAlpha(auto_optimize=True)
    
    # Или ручная настройка с проверкой
    # optimized_alpha = AdaptiveAlpha(1.0, 0.5, 100.0, False, auto_optimize=False)
    
    # В интеграторе добавляем проверку
    def adaptive_system(t, vars):
        n, y, v, h, l = vars
        b = AdaptiveBallistics(n, y, v, h)
        b.optimize_if_needed(t)  # Проверка и оптимизация при необходимости
        return [
            b.delta_polar(t),
            b.delta_trajangle(t),
            b.delta_velocity(t),
            b.delta_altitude(t),
            b.delta_longitude(t)
        ]