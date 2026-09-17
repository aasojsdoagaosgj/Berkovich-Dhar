"""Exact finite coefficients and elementary rational margins for quartic/INFINITE_SIGNS.md."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib, json

BASE = Path(__file__).resolve().parent


def exp_lower(x, n):
    x = Q(x)
    term = total = Q(1)
    for k in range(1, n + 1):
        term *= x / k
        total += term
    return total


def coefficients(M=300):
    a = [1] + [0] * M
    for j in range(1, M + 1):
        if j % 3:
            for _ in range(4):
                for k in range(M, j - 1, -1):
                    a[k] -= a[k - j]
    return a


def run():
    a = coefficients(300)
    signs = [1, -1, 1]
    assert all(v * signs[j % 3] > 0 for j, v in enumerate(a))

    u = Q(157, 900)  # pi/18 > 157/900
    e17 = exp_lower(17, 80)
    checks = {
        'phase_margin': 2 * (u - u**3 / 6) > Q(1, 3),
        'constant_error': 6 * 3**11 + 54 * 3**8 < 10**7,
        'error_ratio_at_300': Q(240000, e17) + Q(6 * 10**11, e17 * e17) < Q(1, 50),
    }
    assert all(checks.values()), checks

    result = {
        'status': 'passed',
        'max_degree': 300,
        'signs_by_residue': signs,
        'coefficient_sha256': hashlib.sha256(str(a).encode()).hexdigest(),
        'rational_checks': checks,
        'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Exact coefficients through degree 300 and elementary scalar inequalities used in quartic/INFINITE_SIGNS.md.',
    }
    (BASE / 'results').mkdir(exist_ok=True)
    (BASE / 'results' / 'infinite_certificate.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    run()
