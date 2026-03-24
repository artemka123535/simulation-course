import random
import numpy as np


class LCG:
    def __init__(self, seed=42):
        self.state = seed
        self.m = 2**32
        self.a = 1664525
        self.c = 1

    def next(self):
        self.state = (self.a * self.state + self.c) % self.m
        return self.state / self.m

N = int(input("Введите размер выборки: "))
lcg = LCG(seed=22)

custom_sample = [lcg.next() for i in range(N)]
builtin_sample = [random.random() for i in range(N)] 

def get_stats(sample):
    mean = np.mean(sample)
    variance = np.var(sample)
    return mean, variance

mean_custom, var_custom = get_stats(custom_sample)
mean_builtin, var_builtin = get_stats(builtin_sample)

theoretical_mean = 0.5
theoretical_var = 1/12

print(f"{'Характеристика':<25} | {'Теоритический':<15} | {'Собственный':<20} | {'Встроенный'}")
print("-" * 80)
print(f"{'Мат. ожидание':<25} | {theoretical_mean:<15.5f} | {mean_custom:<20.5f} | {mean_builtin:.5f}")
print(f"{'Дисперсия':<25} | {theoretical_var:<15.5f} | {var_custom:<20.5f} | {var_builtin:.5f}")
print(f"{'Погрешность мат. ожидания':<25} | {0:<15.5f} | {abs(mean_custom - theoretical_mean):<20.5f} | {abs(mean_builtin - theoretical_mean):.5f}")
print(f"{'Погрешность дисперсии':<25} | {0:<15.5f} | {abs(var_custom - theoretical_var):<20.5f} | {abs(var_builtin - theoretical_var):.5f}")