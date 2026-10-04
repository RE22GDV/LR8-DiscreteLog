"""Реальні послідовні вимірювання без моделювання часу або результатів."""
from argparse import ArgumentParser
from dataclasses import asdict
from datetime import datetime, timezone
from math import isqrt
from pathlib import Path
from random import Random
from statistics import median
import gc
import json
import platform
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'solution'))
from codingame_solution import discrete_log
from discretelog import (BabyTable, bsgs, brute_force, interval_bsgs,
                         pohlig_hellman, pollard_rho, factor, is_prime, primitive_root, Result)


def safe_prime(bits):
    lower = 1 << (bits - 2)
    r = lower + 123
    r += 1 - r % 2
    while not (is_prime(r) and is_prime(2 * r + 1)):
        r += 2
    return 2 * r + 1


def smooth_prime(bits):
    lower, upper = (1 << (bits - 1)), min((1 << bits) - 1, 49999999999)
    candidates = []
    a = 2
    while a < upper:
        b = a
        while b < upper:
            c = b
            while c < upper:
                d = c
                while d < upper:
                    if lower <= d + 1 <= upper:
                        candidates.append(d + 1)
                    d *= 7
                c *= 5
            b *= 3
        a *= 2
    return next(p for p in sorted(set(candidates)) if is_prime(p))


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--seconds', type=int, default=600)
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, default=ROOT / 'docs/results/benchmark.json')
    args = parser.parse_args()
    if args.seconds < 30 or args.repeats < 3:
        parser.error('Потрібно щонайменше 30 секунд бюджету та три повтори')
    started = perf_counter()
    data = {'protocol': {'created_utc': datetime.now(timezone.utc).isoformat(),
                         'python': platform.python_version(), 'system': platform.system(),
                         'processor': 'Intel', 'timer': 'time.perf_counter',
                         'budget_seconds': args.seconds, 'repeats': args.repeats,
                         'seed': 8008, 'gc_during_timing': False,
                         'memory_definition': 'sys.getsizeof(dict) + сума розмірів ключів і значень; не RSS'},
            'groups': [], 'baseline': [], 'structure': [], 'rho': [],
            'tradeoff': [], 'reuse': {}, 'interval': [], 'elgamal': {}}

    def checkpoint():
        if perf_counter() - started >= args.seconds:
            raise TimeoutError('Бюджет експерименту вичерпано')

    def measure(callback, expected=None, repeats=None):
        checkpoint()
        warm = callback()
        if expected is not None:
            assert warm.exponent == expected
        times = []
        result = warm
        for _ in range(repeats or args.repeats):
            checkpoint()
            gc.collect()
            gc.disable()
            try:
                t0 = perf_counter()
                result = callback()
                times.append(perf_counter() - t0)
            finally:
                gc.enable()
            if expected is not None:
                assert result.exponent == expected
        return {'samples_seconds': times, 'median_seconds': median(times),
                'min_seconds': min(times), 'max_seconds': max(times),
                'exponent': result.exponent, 'status': result.status,
                'operations': asdict(result.metrics)}

    groups = {}
    try:
        for bits in [12, 16, 20, 24, 28, 32, 36]:
            for kind in ['smooth', 'safe']:
                checkpoint()
                p = smooth_prime(bits) if kind == 'smooth' else safe_prime(bits)
                g = primitive_root(p)
                assert p.bit_length() == bits and p < 50000000000
                item = {'bits': bits, 'kind': kind, 'p': p, 'g': g,
                        'order': p - 1, 'factors': factor(p - 1)}
                groups[bits, kind] = item
                data['groups'].append(item)
        print('Параметри груп підготовлено', flush=True)

        for bits in [12, 16, 20, 24]:
            group = groups[bits, 'safe']
            p, g, n = group['p'], group['g'], group['order']
            x, h = n - 2, pow(g, n - 2, p)
            data['baseline'].append({'bits': bits, 'p': p, 'x': x,
                                     'brute': measure(lambda: brute_force(g, h, p), x, repeats=3),
                                     'bsgs': measure(lambda: bsgs(g, h, p), x)})
            print(f'Перебір і BSGS: {bits} біт', flush=True)

        for bits in [16, 20, 24, 28, 32, 36]:
            for kind in ['smooth', 'safe']:
                group = groups[bits, kind]
                p, g, n = group['p'], group['g'], group['order']
                x, h = n - 2, pow(g, n - 2, p)
                data['structure'].append({'bits': bits, 'kind': kind, 'p': p,
                                          'largest_prime_factor': max(group['factors']),
                                          'bsgs': measure(lambda: bsgs(g, h, p), x),
                                          'ph': measure(lambda: pohlig_hellman(g, h, p, order=n), x)})
            print(f'Структура порядку: {bits} біт', flush=True)

        for bits in [20, 24, 28, 32, 36]:
            group = groups[bits, 'safe']
            p, g, n = group['p'], pow(group['g'], 2, group['p']), (group['p'] - 1) // 2
            x = n * 3 // 4
            h = pow(g, x, p)
            # Незалежні траєкторії: по одному вимірюванню для кожного з 15 зерен.
            rho_runs = []
            for seed in range(8008, 8023):
                checkpoint()
                gc.disable()
                try:
                    t0 = perf_counter()
                    result = pollard_rho(g, h, p, order=n, seed=seed)
                    elapsed = perf_counter() - t0
                finally:
                    gc.enable()
                assert result.exponent == x
                rho_runs.append({'seed': seed, 'seconds': elapsed, **asdict(result.metrics)})
            table = BabyTable(g, p, order_bound=n)
            data['rho'].append({'bits': bits, 'p': p, 'order': n, 'x': x,
                                'bsgs': measure(lambda: bsgs(g, h, p, order_bound=n), x),
                                'bsgs_table_bytes': table.allocated_bytes(),
                                'runs': rho_runs, 'rho_median_seconds': median(r['seconds'] for r in rho_runs)})
            del table
            print(f'Метод Полларда: {bits} біт, 15 траєкторій', flush=True)

        group = groups[28, 'safe']
        p, g, n = group['p'], group['g'], group['order']
        x, h = n - 2, pow(g, n - 2, p)
        root = isqrt(n - 1) + 1
        for numerator, denominator in [(1, 8), (1, 4), (1, 2), (1, 1), (2, 1), (4, 1), (8, 1)]:
            m = root * numerator // denominator
            table = BabyTable(g, p, baby_steps=m)
            memory = table.allocated_bytes()
            del table
            data['tradeoff'].append({'p': p, 'order': n, 'm': m, 'sqrt_order': root,
                                     'table_bytes': memory,
                                     'timing': measure(lambda: bsgs(g, h, p, baby_steps=m), x)})

        group = groups[32, 'safe']
        p, g, n = group['p'], group['g'], group['order']
        rng = Random(8008)
        exponents = [rng.randrange(1, n) for _ in range(20)]
        targets = [pow(g, x, p) for x in exponents]
        def twenty_fresh():
            result = None
            for target, expected in zip(targets, exponents):
                result = bsgs(g, target, p)
                assert result.exponent == expected
            return result
        table = BabyTable(g, p)
        def twenty_reused():
            result = None
            for target, expected in zip(targets, exponents):
                result = table.solve(target)
                assert result.exponent == expected
            return result
        def build_and_twenty():
            batch = BabyTable(g, p)
            result = None
            for target, expected in zip(targets, exponents):
                result = batch.solve(target)
                assert result.exponent == expected
            return result
        data['reuse'] = {'p': p, 'count': 20, 'exponents': exponents,
                         'fresh': measure(twenty_fresh), 'reuse_with_build': measure(build_and_twenty),
                         'queries_only': measure(twenty_reused),
                         'table_bytes': table.allocated_bytes()}
        del table

        p = next(p for p in range(49999999999, 49999990000, -2) if is_prime(p))
        g, n = primitive_root(p), p - 1
        x, h = n - 2, pow(primitive_root(p), n - 2, p)
        data['upper_bound'] = {'p': p, 'g': g, 'h': h, 'x': x,
                               'bsgs': measure(lambda: bsgs(g, h, p), x),
                               'standalone': measure(lambda: Result(discrete_log(g, h, p), 'found'), x),
                               'table_bytes': BabyTable(g, p).allocated_bytes()}
        for width in [1000, 10000, 1000000, 100000000]:
            lower, upper = x - width + 1, x
            data['interval'].append({'width': width, 'lower': lower, 'upper': upper,
                                      'timing': measure(lambda: interval_bsgs(g, h, p, lower, upper), x)})

        # Навчальне ElGamal: параметри, повідомлення й одноразовий показник власні.
        group = groups[36, 'smooth']
        p, g = group['p'], group['g']
        private, nonce, message = 123456789, 7654321, 2026
        public = pow(g, private, p)
        c1 = pow(g, nonce, p)
        c2 = message * pow(public, nonce, p) % p
        attack = pohlig_hellman(g, public, p)
        recovered = c2 * pow(pow(c1, attack.exponent, p), -1, p) % p
        assert recovered == message
        data['elgamal'] = {'p': p, 'g': g, 'private': private, 'public': public,
                           'nonce': nonce, 'message': message, 'c1': c1, 'c2': c2,
                           'recovered_private': attack.exponent, 'recovered_message': recovered,
                           'attack': measure(lambda: pohlig_hellman(g, public, p), private)}
        data['complete'] = True
    except TimeoutError as exc:
        data['complete'] = False
        data['stop_reason'] = str(exc)
    data['wall_seconds'] = perf_counter() - started
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8', newline='\n')
    print(f"Вимірювання збережено: {args.output}; {data['wall_seconds']:.1f} с", flush=True)
    return 0 if data['complete'] else 3


if __name__ == '__main__':
    raise SystemExit(main())
