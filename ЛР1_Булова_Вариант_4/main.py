"""
Лабораторная работа №1. Вариант 4.
Осипова Анна.

Вариант:
a = 0.25, b = 15, c = 10
eps = 1e-4
критерий для корней: по невязке
x0 = 1
Доп. требование Д3: нормализация аргумента для ln(c).
"""
import math
from typing import Tuple

EPS = 1e-4
N_MAX = 100


def sqrt_newton(a: float, eps: float = EPS, n_max: int = N_MAX) -> Tuple[float, int]:
    """sqrt(a) методом Ньютона. Остановка по невязке |x^2-a| < eps."""
    if not isinstance(a, (int, float)) or not math.isfinite(a):
        raise ValueError("a должен быть конечным числом.")
    if a < 0:
        raise ValueError("a должен быть неотрицательным.")
    if a == 0:
        return 0.0, 0

    x = 1.0  # x0 по варианту
    for n in range(1, n_max + 1):
        x_new = (x + a / x) / 2
        if abs(x_new * x_new - a) < eps:
            return x_new, n
        x = x_new
    raise RuntimeError("Точность не достигнута за N_MAX итераций.")


def cbrt_newton(b: float, eps: float = EPS, n_max: int = N_MAX) -> Tuple[float, int]:
    """∛b методом Ньютона. Остановка по невязке |x^3-b| < eps."""
    if not isinstance(b, (int, float)) or not math.isfinite(b):
        raise ValueError("b должен быть конечным числом.")
    if b == 0:
        return 0.0, 0

    sign = -1.0 if b < 0 else 1.0
    value = abs(b)
    x = 1.0
    for n in range(1, n_max + 1):
        x_new = (2 * x + value / (x * x)) / 3
        result = sign * x_new
        if abs(result ** 3 - b) < eps:
            return result, n
        x = x_new
    raise RuntimeError("Точность не достигнута за N_MAX итераций.")


def ln_series(c: float, eps: float = EPS, n_max: int = 100000) -> Tuple[float, int]:
    """ln(c) рядом 2*sum(y^(2k+1)/(2k+1)), y=(c-1)/(c+1)."""
    if not isinstance(c, (int, float)) or not math.isfinite(c):
        raise ValueError("c должен быть конечным числом.")
    if c <= 0:
        raise ValueError("c должен быть положительным.")

    y = (c - 1) / (c + 1)
    y2 = y * y
    power = y
    total = 0.0

    for k in range(n_max):
        term = 2 * power / (2 * k + 1)
        if abs(term) < eps:
            return total, k
        total += term
        power *= y2
    raise RuntimeError("Точность не достигнута за заданное число членов.")


def ln_series_normalized(c: float, eps: float = EPS, n_max: int = 100000) -> Tuple[float, int, int]:
    """
    Нормализованный вариант Д3.
    c = m*2^k, m в [1, 2), ln(c)=ln(m)+k*ln(2).
    ln(m) и ln(2) вычисляются тем же рядом.
    """
    if not isinstance(c, (int, float)) or not math.isfinite(c):
        raise ValueError("c должен быть конечным числом.")
    if c <= 0:
        raise ValueError("c должен быть положительным.")

    k = math.floor(math.log2(c))
    m = c / (2 ** k)

    # Ряд для ln(m)
    y = (m - 1) / (m + 1)
    y2 = y * y
    power = y
    total_m = 0.0
    count_m = 0
    for j in range(n_max):
        term = 2 * power / (2 * j + 1)
        if abs(term) < eps:
            count_m = j
            break
        total_m += term
        power *= y2
    else:
        raise RuntimeError("Не удалось вычислить ln(m).")

    # Ряд для ln(2), один раз в рамках вычисления.
    y2c = (2 - 1) / (2 + 1)
    power = y2c
    total_2 = 0.0
    count_2 = 0
    for j in range(n_max):
        term = 2 * power / (2 * j + 1)
        if abs(term) < eps:
            count_2 = j
            break
        total_2 += term
        power *= y2c * y2c
    else:
        raise RuntimeError("Не удалось вычислить ln(2).")

    return total_m + k * total_2, count_m + k * count_2, k


def main():
    a, b, c = 0.25, 15.0, 10.0

    s, ns = sqrt_newton(a)
    cb, ncb = cbrt_newton(b)
    l, nl = ln_series(c)
    ln_norm, nnorm, k = ln_series_normalized(c)

    print("Лабораторная работа №1 — Осипова Анна, вариант 4")
    print(f"a={a}, b={b}, c={c}, eps={EPS:g}, критерий: невязка, x0=1")
    print()
    print(f"{'Функция':<18}{'Результат':>18}{'Эталон':>18}{'Погрешность':>18}{'Итерации':>12}")
    print("-" * 84)
    rows = [
        ("sqrt(a)", s, math.sqrt(a), ns),
        ("cbrt(b)", cb, math.copysign(abs(b) ** (1/3), b), ncb),
        ("ln(c)", l, math.log(c), nl),
        ("ln(c), Д3", ln_norm, math.log(c), nnorm),
    ]
    for name, value, ref, count in rows:
        print(f"{name:<18}{value:>18.12f}{ref:>18.12f}{abs(value-ref):>18.3e}{count:>12}")

    print(f"\nД3: c = m * 2^k -> m = {c/(2**k):.6f}, k = {k}")
    print(f"Ненормированный ln(c): {nl} итераций")
    print(f"Нормированный ln(c):   {nnorm} итераций")


if __name__ == "__main__":
    main()
