# Six eventual Berkovich–Dhar cases

This directory contains the main proof drafts and computational checks for

```math
(p,s)=(3,5),(3,6),(3,7),(3,8),(5,2),(5,3).
```

The claims are **eventual**: for each case, the stated coefficient sign pattern and transition description hold for all sufficiently large $n$. No explicit threshold for the six-case theorem is claimed here, and the original all-$n$ conjectures are not claimed proved.

## Proof files

- `MOD3_RESULT.md`: modulus 3, powers 5–8.
- `MOD5_RESULT.md`: modulus 5, powers 2–3 and the global range connection.
- `MOD5_COMPACT.md`: compact saddle analysis and transition uniqueness.
- `MOD5_ENDPOINT.md`: growing-saddle analysis, including vanishing cube residues.
- `MOD5_FIRST_TAIL.md`: first nonzero tail of the modulus-5 cube.
- `INFINITE_SIGNS.md`: infinite-product sign input.
- `AUDIT.md`: original internal completion checklist (not independent peer review).

## Computations

`verify.py` checks saved reports and their source hashes. With `--replay` it reruns the finite diagnostics and the retained modulus-5 interval certificates.

```bash
python -B berkovich-dhar-six/verify.py
python -B berkovich-dhar-six/verify.py --replay
```

The required Python package is pinned in the repository-root `requirements.txt`.

The JSON files under `results/` are committed so that the exact outputs used in the draft can be inspected. They are supporting certificates/diagnostics, not a formal proof checker for the analytic arguments.
