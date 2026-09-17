# Reproducibility and computation scope

The analytical arguments in the Markdown proof files are ordinary mathematical proofs. The accompanying Python programs serve three narrower purposes: exact finite coefficient checks, interval/rational verification of scalar inequalities, and numerical diagnostics of asymptotic expansions.

## One-command checks

From the repository root:

```bash
python -m pip install -r requirements.txt
python -B berkovich-dhar-six/verify.py
```

This checks that the saved JSON reports match the committed scripts and dependencies and that the expected six cases are present.

To regenerate the computational reports:

```bash
python -B berkovich-dhar-six/verify.py --replay
```

Do not use `python -O`.

## Main scripts

| Script | Scope |
|---|---|
| `berkovich-dhar-six/certify_infinite.py` | Exact coefficients through degree 2000 for all six associated infinite products, plus rational checks used with Wang's explicit remainder bound. |
| `berkovich-dhar-six/certify_first_tail.py` | Exact polynomial identities and truncated diagnostics for the positive-kernel argument in the modulus-5 cube. The infinite coefficientwise inequalities come from the written positive expansions, not the truncation. |
| `berkovich-dhar-six/finite_diagnostic.py` | Every lower-half coefficient for each of the six finite polynomials for $n=1,\dots,40$, using exact integers. Diagnostic only. |
| `berkovich-dhar-six/saddle_diagnostic.py` | Four coefficients near each of the eight transitions at $n=20,40,80$, compared to the two-term saddle correction at high precision. Diagnostic only. |
| `berkovich-dhar-six/endpoint_diagnostic.py` | Rational scalar checks and high-precision direct comparisons for the modulus-5 endpoint formulas. |
| `third-borwein-mod5/verification/certify_scalars.py` | Re-runs the retained interval certificates supporting the modulus-5 localisation dependency. |
| `quartic/certify_infinite.py` | Exact degree-0-through-300 certificate and scalar inequalities for the modulus-3 quartic infinite-product dependency. |
| `quartic/certify_endpoint_constants.py` | Scalar margin checks for the modulus-3 quartic endpoint dependency. |

## Saved reports

`berkovich-dhar-six/results/` contains the outputs used by `verify.py`, including:

- `infinite_certificate.json`
- `first_tail_certificate.json`
- `finite_diagnostic.json`
- `saddle_diagnostic.json`
- `endpoint_diagnostic.json`
- `profiles.json`
- `dependencies/*.json`

The largest files are committed deliberately so the interval-certified dependency checks are inspectable without regenerating them first.

## Important limitation

Passing every script does **not** mechanically prove the large-$n$ theorem. In particular, uniformity of the Euler–Maclaurin/saddle expansions, the complementary-arc analysis, and the logical connection among coefficient ranges remain mathematical arguments in the proof text. The second-pass AI audit specifically reviewed those interfaces; see `audit/AI_AUDIT_20260918.md`.
