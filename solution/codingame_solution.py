"""Самодостатній розв'язок Discrete Log Problem для CodinGame."""
from math import isqrt
import sys


def discrete_log(g, h, q):
    # Q просте; показники достатньо перевірити в межах 0..Q-2.
    n = q - 1
    m = isqrt(n - 1) + 1
    baby = {}
    value = 1
    for j in range(m):
        baby.setdefault(value, j)  # Зберігаємо найменший показник.
        value = value * g % q

    factor = pow(value, -1, q)  # G^(-m) mod Q.
    gamma = h
    for i in range((n + m - 1) // m):
        j = baby.get(gamma)
        if j is not None:
            x = i * m + j
            if x < n:
                return x
        gamma = gamma * factor % q
    return None


def main():
    g, h, q = map(int, sys.stdin.read().split())
    x = discrete_log(g, h, q)
    if x is None:
        raise ValueError("У заданій підгрупі розв'язку немає")
    print(x)


if __name__ == "__main__":
    main()
