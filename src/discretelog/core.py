"""Пошук найменшого невід'ємного дискретного логарифма в простому полі."""
from dataclasses import dataclass, field
from math import isqrt
from random import Random
import sys

from .number_theory import crt, factor, is_prime, multiplicative_order, validate


@dataclass
class Metrics:
    baby_multiplications: int = 0
    giant_multiplications: int = 0
    lookups: int = 0
    table_entries: int = 0
    brute_multiplications: int = 0
    rho_transitions: int = 0
    rho_restarts: int = 0
    ph_digits: int = 0


@dataclass
class Result:
    exponent: int | None
    status: str
    metrics: Metrics = field(default_factory=Metrics)


def brute_force(g, h, p, *, max_steps=None):
    g, h, p = validate(g, h, p)
    if max_steps is not None and max_steps < 0:
        raise ValueError("Бюджет не може бути від'ємним")
    metrics = Metrics()
    value = 1
    limit = p - 1 if max_steps is None else min(max_steps, p - 1)
    for x in range(limit):
        if value == h:
            return Result(x, "found", metrics)
        value = value * g % p
        metrics.brute_multiplications += 1
        if value == 1:
            return Result(None, "no_solution", metrics)
    return Result(None, "no_solution" if limit == p - 1 else "budget_exhausted", metrics)


class BabyTable:
    """Таблиця для повторних запитів з незмінними G, P і межею порядку."""
    def __init__(self, g, p, *, order_bound=None, baby_steps=None):
        self.g, _, self.p = validate(g, 1, p)
        self.n = p - 1 if order_bound is None else order_bound
        if self.n < 1 or pow(self.g, self.n, p) != 1:
            raise ValueError("Межа має бути додатним кратним порядку G")
        self.m = isqrt(self.n - 1) + 1 if baby_steps is None else baby_steps
        if not isinstance(self.m, int) or isinstance(self.m, bool) or not 1 <= self.m <= self.n:
            raise ValueError("Розмір таблиці має належати 1..order_bound")
        self.baby = {}
        value = 1
        for j in range(self.m):
            self.baby.setdefault(value, j)
            value = value * self.g % p
        self.factor = pow(value, -1, p)

    def allocated_bytes(self):
        # Облік Python-об'єктів; це не RSS і не пікова пам'ять процесу.
        return sys.getsizeof(self.baby) + sum(sys.getsizeof(k) + sys.getsizeof(v) for k, v in self.baby.items())

    def solve(self, h):
        if not isinstance(h, int) or isinstance(h, bool) or h % self.p == 0:
            raise ValueError("Потрібен ненульовий цілий елемент")
        h %= self.p
        metrics = Metrics(table_entries=len(self.baby))
        gamma = h
        for i in range((self.n + self.m - 1) // self.m):
            metrics.lookups += 1
            j = self.baby.get(gamma)
            if j is not None:
                x = i * self.m + j
                if x < self.n and pow(self.g, x, self.p) == h:
                    return Result(x, "found", metrics)
            gamma = gamma * self.factor % self.p
            metrics.giant_multiplications += 1
        return Result(None, "no_solution", metrics)


def bsgs(g, h, p, *, order_bound=None, baby_steps=None):
    g, h, p = validate(g, h, p)
    table = BabyTable(g, p, order_bound=order_bound, baby_steps=baby_steps)
    result = table.solve(h)
    result.metrics.baby_multiplications = table.m
    return result


def interval_bsgs(g, h, p, lower, upper):
    g, h, p = validate(g, h, p)
    if not 0 <= lower <= upper:
        raise ValueError("Потрібен непорожній невід'ємний інтервал")
    width = upper - lower + 1
    m = isqrt(width - 1) + 1
    baby, value = {}, 1
    for j in range(m):
        baby.setdefault(value, j)
        value = value * g % p
    gamma = h * pow(pow(g, lower, p), -1, p) % p
    step = pow(value, -1, p)
    metrics = Metrics(baby_multiplications=m, table_entries=len(baby))
    for i in range((width + m - 1) // m):
        metrics.lookups += 1
        j = baby.get(gamma)
        if j is not None:
            offset = i * m + j
            if offset < width:
                return Result(lower + offset, "found", metrics)
        gamma = gamma * step % p
        metrics.giant_multiplications += 1
    return Result(None, "no_solution", metrics)


def pohlig_hellman(g, h, p, *, order=None):
    g, h, p = validate(g, h, p)
    n = multiplicative_order(g, p) if order is None else order
    factors = factor(n)
    if n < 1 or pow(g, n, p) != 1 or any(pow(g, n // prime, p) == 1 for prime in factors):
        raise ValueError("Потрібен точний порядок G")
    metrics = Metrics()
    if pow(h, n, p) != 1:
        return Result(None, "no_solution", metrics)
    if n == 1:
        return Result(0 if h == 1 else None, "found" if h == 1 else "no_solution", metrics)
    residues, moduli = [], []
    inverse_g = pow(g, -1, p)
    for prime, power in factors.items():
        digit_base = pow(g, n // prime, p)
        table = BabyTable(digit_base, p, order_bound=prime)
        metrics.baby_multiplications += table.m
        metrics.table_entries = max(metrics.table_entries, len(table.baby))
        residue, prime_power = 0, 1
        for _ in range(power):
            target = pow(h * pow(inverse_g, residue, p) % p, n // (prime_power * prime), p)
            digit = table.solve(target)
            metrics.ph_digits += 1
            metrics.giant_multiplications += digit.metrics.giant_multiplications
            metrics.lookups += digit.metrics.lookups
            if digit.exponent is None:
                return Result(None, "no_solution", metrics)
            residue += digit.exponent * prime_power
            prime_power *= prime
        residues.append(residue)
        moduli.append(prime_power)
    x = crt(residues, moduli)
    if pow(g, x, p) != h:
        raise AssertionError("Перевірка результату CRT не пройшла")
    return Result(x, "found", metrics)


def pollard_rho(g, h, p, *, order, seed=8008, max_restarts=16, max_iterations=None):
    """16 розгалужень та пошук циклу Флойда; лише простий порядок підгрупи."""
    g, h, p = validate(g, h, p)
    if not is_prime(order) or g == 1 or pow(g, order, p) != 1:
        raise ValueError("Метод стенда потребує простого порядку G")
    metrics = Metrics()
    if pow(h, order, p) != 1:
        return Result(None, "no_solution", metrics)
    if h == 1:
        return Result(0, "found", metrics)
    random = Random(seed)
    limit = max_iterations if max_iterations is not None else 12 * (isqrt(order) + 1)
    for attempt in range(max_restarts):
        metrics.rho_restarts = attempt
        increments = []
        for _ in range(16):
            a, b = random.randrange(order), random.randrange(order)
            increments.append((pow(g, a, p) * pow(h, b, p) % p, a, b))

        def step(state):
            value, a, b = state
            multiplier, da, db = increments[value % 16]
            metrics.rho_transitions += 1
            return value * multiplier % p, (a + da) % order, (b + db) % order

        a, b = random.randrange(order), random.randrange(order)
        slow = fast = (pow(g, a, p) * pow(h, b, p) % p, a, b)
        for _ in range(limit):
            slow, fast = step(slow), step(step(fast))
            if slow[0] == fast[0]:
                denominator = (fast[2] - slow[2]) % order
                if denominator:
                    x = (slow[1] - fast[1]) * pow(denominator, -1, order) % order
                    if pow(g, x, p) == h:
                        return Result(x, "found", metrics)
                break
    return Result(None, "budget_exhausted", metrics)
