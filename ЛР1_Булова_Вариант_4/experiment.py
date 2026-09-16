"""Дополнительный эксперимент Д3: сравнение ln(c) до и после нормализации."""
import math
from main import ln_series, ln_series_normalized

if __name__ == "__main__":
    c = 10.0
    eps = 1e-4
    raw, n_raw = ln_series(c, eps)
    norm, n_norm, k = ln_series_normalized(c, eps)
    ref = math.log(c)

    print("Д3 — сравнение способов вычисления ln(c)")
    print(f"c={c}, eps={eps}")
    print(f"Ненормированный: {raw:.12f}, итераций={n_raw}, ошибка={abs(raw-ref):.3e}")
    print(f"Нормированный:   {norm:.12f}, итераций={n_norm}, ошибка={abs(norm-ref):.3e}")
    print(f"Нормализация: c = m*2^k, k={k}, m={c/(2**k):.6f}")
