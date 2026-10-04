"""Побудова PNG та векторних PDF виключно зі збережених вимірювань."""
from pathlib import Path
import json
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'docs/results/benchmark.json'
OUT = ROOT / 'docs/figures'
BLUE, ORANGE, GREEN, GRAY = '#27639b', '#b35b22', '#278066', '#555555'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.grid': True, 'grid.alpha': .22, 'pdf.fonttype': 42})


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'pdf').mkdir(exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT / f'{name}.png', dpi=180)
    fig.savefig(OUT / 'pdf' / f'{name}.png', dpi=180)
    fig.savefig(OUT / 'pdf' / f'{name}.pdf')
    plt.close(fig)


def main():
    d = json.loads(DATA.read_text(encoding='utf-8'))
    if not d['complete']:
        raise ValueError('Потрібна завершена серія експериментів')
    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    for key, label, color in [('brute', 'Повний перебір', ORANGE), ('bsgs', 'BSGS', BLUE)]:
        rows = d['baseline']
        x = [r['bits'] for r in rows]
        ax.plot(x, [r[key]['median_seconds'] * 1000 for r in rows], 'o-', color=color, label=label)
        ax.fill_between(x, [r[key]['min_seconds'] * 1000 for r in rows],
                        [r[key]['max_seconds'] * 1000 for r in rows], color=color, alpha=.15)
    ax.set(xlabel='Довжина модуля, біт', ylabel='Час, мс (логарифмічна шкала)', yscale='log')
    ax.set_xticks([12, 16, 20, 24]); ax.legend()
    save(fig, '01_brute_bsgs')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    for kind, label, color in [('smooth', 'Гладкий порядок', GREEN), ('safe', 'Порядок 2r, r просте', ORANGE)]:
        rows = [r for r in d['structure'] if r['kind'] == kind]
        x = [r['bits'] for r in rows]
        ax.plot(x, [r['ph']['median_seconds'] * 1000 for r in rows], 'o-', color=color, label=f'Поліг—Геллман: {label}')
        ax.plot(x, [r['bsgs']['median_seconds'] * 1000 for r in rows], '--', color=color, alpha=.7, label=f'BSGS: {label}')
    ax.set(xlabel='Довжина модуля, біт', ylabel='Час, мс (логарифмічна шкала)', yscale='log')
    ax.set_xticks([16, 20, 24, 28, 32, 36]); ax.legend(fontsize=9)
    save(fig, '02_order_structure')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    for kind, label, color in [('smooth', 'Гладкий порядок', GREEN), ('safe', 'Порядок 2r', ORANGE)]:
        rows = [r for r in d['structure'] if r['kind'] == kind]
        x = [r['largest_prime_factor'] for r in rows]
        y = [r['ph']['median_seconds'] * 1000 for r in rows]
        ax.scatter(x, y, color=color, s=42, label=label)
        for r, xx, yy in zip(rows, x, y):
            if kind == 'smooth' and r['bits'] not in (16, 36):
                continue
            ax.annotate(f"{r['bits']}", (xx, yy), xytext=(5, 3), textcoords='offset points', fontsize=8)
    ax.set(xlabel='Найбільший простий множник порядку', ylabel='Поліг—Геллман, мс', xscale='log', yscale='log')
    ax.legend()
    save(fig, '03_largest_factor')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    rows = d['rho']; x = [r['bits'] for r in rows]
    ax.plot(x, [r['bsgs']['median_seconds'] * 1000 for r in rows], 'o-', color=BLUE, label='BSGS')
    ax.plot(x, [r['rho_median_seconds'] * 1000 for r in rows], 'o-', color=ORANGE, label='Поллард ρ, медіана 15 траєкторій')
    ax.fill_between(x, [np.percentile([t['seconds'] * 1000 for t in r['runs']], 25) for r in rows],
                    [np.percentile([t['seconds'] * 1000 for t in r['runs']], 75) for r in rows], color=ORANGE, alpha=.18, label='Міжквартильний діапазон')
    ax.set(xlabel='Довжина модуля, біт', ylabel='Час, мс (логарифмічна шкала)', yscale='log')
    ax.set_xticks(x); ax.legend(fontsize=9)
    save(fig, '04_pollard_bsgs')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    rows = d['tradeoff']; x = [r['m'] / r['sqrt_order'] for r in rows]
    time_line, = ax.plot(x, [r['timing']['median_seconds'] * 1000 for r in rows], 'o-', color=BLUE)
    ax.set(xlabel='Розмір таблиці / ⌈√n⌉', ylabel='Час BSGS, мс', xscale='log', xticks=x)
    ax.set_xticklabels(['1/8', '1/4', '1/2', '1', '2', '4', '8'])
    ax2 = ax.twinx(); ax2.spines['right'].set_visible(True); ax2.grid(False)
    memory_line, = ax2.plot(x, [r['table_bytes'] / 2**20 for r in rows], 's--', color=GRAY)
    ax2.set_ylabel('Облік об’єктів таблиці, MiB', color=GRAY)
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.legend([time_line, memory_line], ['Час', 'Облік об’єктів таблиці'], loc='upper center', fontsize=10)
    save(fig, '05_memory_tradeoff')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    keys = ['fresh', 'reuse_with_build', 'queries_only']
    values = [d['reuse'][key]['median_seconds'] * 1000 for key in keys]
    bars = ax.bar(['20 нових таблиць', 'Одна таблиця + 20 запитів', 'Лише 20 запитів'], values, color=[ORANGE, BLUE, GREEN], width=.65)
    for bar, value in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, value + 3, f'{value:.1f}', ha='center')
    ax.set(ylabel='Сумарний час 20 запитів, мс', ylim=(0, max(values) * 1.2))
    ax.tick_params(axis='x', labelsize=9)
    save(fig, '06_precomputation')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    rows = d['interval']
    ax.plot([r['width'] for r in rows], [r['timing']['median_seconds'] * 1000 for r in rows], 'o-', color=GREEN, label='BSGS у відомому інтервалі')
    ax.axhline(d['upper_bound']['bsgs']['median_seconds'] * 1000, color=BLUE, linestyle='--', label='BSGS у всій групі')
    ax.set(xlabel='Кількість можливих показників у відомому інтервалі', ylabel='Час, мс', xscale='log', yscale='log')
    ax.legend(fontsize=9)
    save(fig, '07_known_interval')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    rows = d['baseline']; x = [r['p'] - 1 for r in rows]
    ax.plot(x, [r['brute']['operations']['brute_multiplications'] for r in rows], 'o-', color=ORANGE, label='Множення повного перебору')
    ax.plot(x, [r['bsgs']['operations']['baby_multiplications'] + r['bsgs']['operations']['giant_multiplications'] for r in rows], 'o-', color=BLUE, label='Множення малих і великих кроків')
    ax.set(xlabel='Порядок групи n', ylabel='Кількість множень (без pow та перевірок)', xscale='log', yscale='log')
    ax.legend(fontsize=9)
    save(fig, '08_operation_counts')

    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    runs = d['rho'][-1]['runs']
    ax.bar([str(r['seed']) for r in runs], [r['seconds'] * 1000 for r in runs], color=ORANGE)
    ax.axhline(statistics.median(r['seconds'] * 1000 for r in runs), color=GRAY, linestyle='--', label='Медіана')
    ax.set(xlabel='Зерно випадкової траєкторії', ylabel='Поллард ρ у 36-бітному модулі, мс')
    ax.tick_params(axis='x', rotation=45, labelsize=8); ax.legend()
    save(fig, '09_pollard_variation')
    print('Побудовано 9 графіків у PNG та векторному PDF')


if __name__ == '__main__':
    main()
