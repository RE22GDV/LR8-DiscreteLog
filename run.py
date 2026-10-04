"""Український командний інтерфейс дослідницького стенда."""
from pathlib import Path
import argparse
from dataclasses import asdict
import json
import sys
from time import perf_counter

sys.path.insert(0, str(Path(__file__).resolve().parent/'src'))
from discretelog import bsgs, brute_force, multiplicative_order, pohlig_hellman, pollard_rho


def main():
    parser = argparse.ArgumentParser(description='Відновлення показника G^X mod P = H')
    sub = parser.add_subparsers(dest='command', required=True)
    solve = sub.add_parser('solve', help='Розв’язати дискретний логарифм')
    for name in ('g','h','p'): solve.add_argument(name, type=int)
    solve.add_argument('--method', choices=['bsgs','brute','ph','rho'], default='bsgs')
    solve.add_argument('--max-steps', type=int, default=1_000_000)
    solve.add_argument('--baby-steps', type=int)
    solve.add_argument('--json', action='store_true')
    sub.add_parser('selftest', help='Перевірити контрольні приклади')
    args = parser.parse_args()
    if args.command == 'selftest':
        for g,h,p,x in [(654,4547,11087,114),(2,8,13,3),(4,2,7,2),(2,1,13,0)]:
            assert bsgs(g,h,p).exponent == x
            assert pohlig_hellman(g,h,p).exponent == x
        print('Контрольні приклади пройдено: 4 з 4.')
        return 0
    started = perf_counter()
    try:
        if args.method == 'brute': result = brute_force(args.g,args.h,args.p,max_steps=args.max_steps)
        elif args.method == 'ph': result = pohlig_hellman(args.g,args.h,args.p)
        elif args.method == 'rho':
            result = pollard_rho(args.g,args.h,args.p,order=multiplicative_order(args.g,args.p))
        else: result = bsgs(args.g,args.h,args.p,baby_steps=args.baby_steps)
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    elapsed = perf_counter()-started
    if args.json:
        print(json.dumps({**asdict(result),'seconds':elapsed},ensure_ascii=False))
    else:
        print(f'Показник: {result.exponent}\nСтан: {result.status}\nЧас: {elapsed:.6f} с')
    return 0 if result.status == 'found' else 3 if result.status == 'budget_exhausted' else 1


if __name__ == '__main__':
    raise SystemExit(main())
