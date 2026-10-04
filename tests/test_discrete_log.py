"""Незалежна перевірка за повною таблицею степенів та синтетичними групами."""
from random import Random
import subprocess
import sys
from pathlib import Path

import pytest

from codingame_solution import discrete_log
from discretelog import (BabyTable, bsgs, brute_force, interval_bsgs,
                         pohlig_hellman, pollard_rho, is_prime,
                         factor, multiplicative_order, primitive_root, crt)

PRIMES = [p for p in range(2, 80) if is_prime(p)]


def oracle(g, p):
    answers, value = {}, 1
    for x in range(p - 1):
        answers.setdefault(value, x)
        value = value * g % p
    return answers


@pytest.mark.parametrize('p', PRIMES)
def test_all_small_groups_against_independent_oracle(p):
    for g in range(1, p):
        answers = oracle(g, p)
        for h in range(1, p):
            expected = answers.get(h)
            for solve in (discrete_log,):
                assert solve(g, h, p) == expected
            for solve in (bsgs, pohlig_hellman, brute_force):
                result = solve(g, h, p)
                assert result.exponent == expected
                assert result.status == ('found' if expected is not None else 'no_solution')


@pytest.mark.parametrize('p', [97, 193, 257, 769, 12289, 65537, 11087, 1000000007])
def test_orders_prime_powers_and_random_exponents(p):
    rng = Random(8008 + p)
    generator = primitive_root(p)
    for _ in range(12):
        g = pow(generator, rng.randrange(1, p - 1), p)
        order = multiplicative_order(g, p)
        x = rng.randrange(0, p - 1)
        h = pow(g, x, p)
        assert pohlig_hellman(g, h, p).exponent == x % order
        assert bsgs(g, h, p, order_bound=order).exponent == x % order


@pytest.mark.parametrize('m', [1, 2, 3, 7, 17, 96])
def test_baby_table_size_and_reuse(m):
    table = BabyTable(5, 97, baby_steps=m)
    for x in range(96):
        assert table.solve(pow(5, x, 97)).exponent == x
    assert table.allocated_bytes() > 0


@pytest.mark.parametrize('lower,upper', [(0, 0), (0, 5), (8, 25), (95, 130), (100, 105)])
def test_interval_returns_smallest_exponent_in_interval(lower, upper):
    for h in range(1, 97):
        expected = next((x for x in range(lower, upper + 1) if pow(5, x, 97) == h), None)
        assert interval_bsgs(5, h, 97, lower, upper).exponent == expected


@pytest.mark.parametrize('p', [7, 23, 47, 167, 1019, 1000000007])
def test_pollard_prime_subgroup_with_restarts(p):
    order = (p - 1) // 2
    g = pow(primitive_root(p), 2, p)
    for seed in range(4):
        x = 1 + seed * (order - 1) // 4
        result = pollard_rho(g, pow(g, x, p), p, order=order, seed=seed)
        assert result.status == 'found'
        assert result.exponent == x


def test_budget_is_distinct_from_no_solution():
    assert brute_force(2, 8, 13, max_steps=2).status == 'budget_exhausted'
    assert brute_force(4, 3, 7).status == 'no_solution'
    assert pollard_rho(4, 2, 7, order=3, max_iterations=0).status == 'budget_exhausted'
    assert pollard_rho(4, 3, 7, order=3).status == 'no_solution'


@pytest.mark.parametrize('g,h,p', [(0, 2, 7), (2, 0, 7), (2, 3, 9), (2, 3, 1), (True, 3, 7)])
def test_invalid_group_inputs(g, h, p):
    with pytest.raises(ValueError):
        bsgs(g, h, p)


def test_invalid_algorithm_parameters():
    with pytest.raises(ValueError): BabyTable(5, 97, baby_steps=0)
    with pytest.raises(ValueError): BabyTable(5, 97, order_bound=95)
    with pytest.raises(ValueError): pohlig_hellman(4, 2, 7, order=6)
    with pytest.raises(ValueError): pollard_rho(5, 8, 97, order=96)
    with pytest.raises(ValueError): interval_bsgs(5, 8, 97, 10, 9)
    with pytest.raises(ValueError): brute_force(5, 8, 97, max_steps=-1)


def test_number_theory():
    assert factor(1) == {}
    assert factor(2**8 * 3**3 * 17) == {2: 8, 3: 3, 17: 1}
    assert not is_prime(341550071728321)
    assert is_prime(49999999967)
    assert crt([2, 3, 2], [3, 5, 7]) == 23
    with pytest.raises(ValueError): crt([1, 1], [2, 4])


def test_platform_stdin_stdout_and_example():
    script = Path(__file__).resolve().parents[1] / 'solution/codingame_solution.py'
    process = subprocess.run([sys.executable, str(script)], input='654 4547 11087\n',
                             text=True, capture_output=True, check=True)
    assert process.stdout == '114\n'
    assert process.stderr == ''


def test_near_platform_upper_bound():
    p = 49999999967
    g = primitive_root(p)
    x = p - 2
    h = pow(g, x, p)
    assert discrete_log(g, h, p) == x
    assert bsgs(g, h, p).exponent == x
    assert pohlig_hellman(g, h, p).exponent == x
