"""Цілочислова арифметика для навчальних простих модулів до 2**64."""
from math import prod


def is_prime(n):
    if not isinstance(n, int) or isinstance(n, bool):
        return False
    if n < 2:
        return False
    if n >= 2**64:
        raise ValueError("Детермінована перевірка обмежена 64 бітами")
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        a %= n
        if a == 0:
            continue
        value = pow(a, d, n)
        if value in (1, n - 1):
            continue
        for _ in range(s - 1):
            value = value * value % n
            if value == n - 1:
                break
        else:
            return False
    return True


def factor(n):
    if not isinstance(n, int) or isinstance(n, bool) or n < 1:
        raise ValueError("Потрібне додатне ціле число")
    factors = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def validate(g, h, p):
    if any(not isinstance(x, int) or isinstance(x, bool) for x in (g, h, p)):
        raise ValueError("Параметри мають бути цілими числами")
    if not is_prime(p):
        raise ValueError("Потрібен простий модуль")
    g, h = g % p, h % p
    if g == 0 or h == 0:
        raise ValueError("Потрібні ненульові елементи мультиплікативної групи")
    return g, h, p


def multiplicative_order(g, p):
    g, _, p = validate(g, 1, p)
    n = p - 1
    for prime in factor(n):
        while n % prime == 0 and pow(g, n // prime, p) == 1:
            n //= prime
    return n


def primitive_root(p):
    if not is_prime(p):
        raise ValueError("Потрібен простий модуль")
    primes = list(factor(p - 1))
    for g in range(1, p):
        if all(pow(g, (p - 1) // prime, p) != 1 for prime in primes):
            return g
    raise AssertionError("У простого модуля має бути первісний корінь")


def crt(residues, moduli):
    if len(residues) != len(moduli) or not moduli:
        raise ValueError("Потрібні узгоджені непорожні набори")
    modulus = prod(moduli)
    result = 0
    for residue, m in zip(residues, moduli):
        if m <= 1:
            raise ValueError("Модулі мають бути більшими за одиницю")
        cofactor = modulus // m
        result += residue * cofactor * pow(cofactor, -1, m)
    return result % modulus
