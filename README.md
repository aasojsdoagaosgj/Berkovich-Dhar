# Eventual sign patterns for six Berkovich–Dhar Borwein-type cases

This repository contains a proof draft for the **eventual (large-$n$) form** of six cases from Berkovich and Dhar's *New Borwein-Type Conjectures*:

```math
(p,s)=(3,5),(3,6),(3,7),(3,8),(5,2),(5,3).
```

The result proved here is **not** the full all-$n$ conjecture.  The claim is that, for each of the six cases, there exists an integer threshold $N_0$ such that the conjectured coefficient-sign pattern holds for every $n\ge N_0$, together with convergence (and a linear correction) for the sign-transition location.

> **Status.** Research draft, 18 September 2026. The derivation, exposition, code, and computational checks were produced with OpenAI models. A separate second-pass AI audit found no obstruction to the stated eventual theorems, but this is **not human peer review or formal verification**. See [`audit/AI_AUDIT_20260918.md`](audit/AI_AUDIT_20260918.md).

## Main results

For modulus $3$, powers $s=5,6,7,8$, the proof establishes eventual positivity in residue $0$, a single sign transition in residue $2$, and the remaining signs by reciprocity.  If $k_{n,s}$ denotes the transition index, then

```math
k_{n,s}=\alpha_s n^2+\beta_s n+O(1).
```

For modulus $5$, powers $s=2,3$, residues $3$ and $4$ each have an eventual single transition (with the weak inequalities required by the original conjecture), and

```math
k_{n,s,a}=\alpha_{s,a}n^2+\beta_{s,a}n+O(1).
```

The limiting constants obtained from the exact phase equations are:

| modulus | power | transition residue | $\alpha$ |
|---:|---:|---:|---:|
| 3 | 5 | 2 | 1.269283684929835... |
| 3 | 6 | 2 | 1.777289033919982... |
| 3 | 7 | 2 | 2.281420279605300... |
| 3 | 8 | 2 | 2.783868331238400... |
| 5 | 2 | 3 | 1.135460959422580... |
| 5 | 2 | 4 | 1.005280823590582... |
| 5 | 3 | 3 | 2.136899919694871... |
| 5 | 3 | 4 | 1.984296232170181... |

Several of these differ from the decimal values printed in the original conjecture paper.  In this repository the limits are defined by exact monotone phase equations; the decimals above are high-precision evaluations, not interval certificates.

The detailed statements are in:

- [`berkovich-dhar-six/MOD3_RESULT.md`](berkovich-dhar-six/MOD3_RESULT.md) — modulus 3, powers 5 through 8.
- [`berkovich-dhar-six/MOD5_RESULT.md`](berkovich-dhar-six/MOD5_RESULT.md) — modulus 5, powers 2 and 3.

## Proof architecture

The proof splits the coefficient range into initial, growing-saddle/endpoint, and compact-saddle regimes.

- [`berkovich-dhar-six/INFINITE_SIGNS.md`](berkovich-dhar-six/INFINITE_SIGNS.md): signs of the associated infinite products, using Liuquan Wang's explicit asymptotic formula plus exact finite certificates.
- [`berkovich-dhar-six/MOD3_RESULT.md`](berkovich-dhar-six/MOD3_RESULT.md): extends the compact and endpoint machinery from the quartic modulus-3 case to powers 5–8.
- [`berkovich-dhar-six/MOD5_COMPACT.md`](berkovich-dhar-six/MOD5_COMPACT.md): four-root saddle amplitude, exact cosine factorisation, unique phase zeros, $O(n^{-3})$ coefficient expansion, and adjacent-coefficient comparison.
- [`berkovich-dhar-six/MOD5_ENDPOINT.md`](berkovich-dhar-six/MOD5_ENDPOINT.md): uniform growing-saddle estimates; in the cube's vanishing residue classes the subtraction is made on the whole coefficient circle so the small tail factor survives every error term.
- [`berkovich-dhar-six/MOD5_FIRST_TAIL.md`](berkovich-dhar-six/MOD5_FIRST_TAIL.md): coefficientwise-positive Rogers–Ramanujan kernels for the first nonzero tail of the modulus-5 cube.
- [`quartic/`](quartic/): modulus-3 quartic compact/endpoint dependencies used by the powers 5–8 argument.
- [`third-borwein-mod5/`](third-borwein-mod5/): the modulus-5 analytic localisation machinery and its interval certificates, from the earlier third Borwein work and retained here as a dependency.  The manuscript there is that dependency, not the statement of the six theorems in this repository.

A key technical point near a transition is that the proof compares two separate value expansions

```math
\frac{c_j}{B_{n,j}}=\Psi(\tau)+\frac{A_1(\tau)}n+\frac{A_2(\tau)}{n^2}+O(n^{-3}),
```

so that for consecutive indices in one residue class

```math
\Delta\!\left(\frac{c_j}{B_{n,j}}\right)
=-\frac{\Psi'(\tau)}{n^2H''(\tau)}+O(n^{-3}).
```

No differentiation of an unspecified $O(n^{-3})$ remainder is used.

## What is and is not claimed

**Claimed in the draft:** existence of a threshold after which all coefficients in each of the six cases have the stated sign behaviour; a unique transition in the relevant residue classes; exact equations for the limiting transition constants; and a refinement $\alpha n^2+\beta n+O(1)$.

**Not claimed:** the original all-positive-$n$ conjectures; a numerical common threshold for these six theorems; formal verification; human peer review; or a mathematically exhaustive priority/novelty determination.

The original conjectures are for every positive integer $n$. Therefore this repository should be cited as an **eventual / large-$n$ result**, not as a complete proof of Berkovich–Dhar Conjectures 2.1 and 2.2.

## Reproduction

Python 3.10 or later is recommended.

```bash
python -m venv .venv
# Windows: .venv\\Scripts\\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt
```

Check the saved computational records and provenance:

```bash
python -B berkovich-dhar-six/verify.py
```

Re-run the finite diagnostics and interval certificates:

```bash
python -B berkovich-dhar-six/verify.py --replay
```

Do **not** run Python with `-O`; the verification scripts intentionally use assertions.  These programs check identities, coefficients, interval inequalities, saved-result provenance, and numerical diagnostics.  They do **not** mechanically verify the analytic proof text.

See [`REPRODUCIBILITY.md`](REPRODUCIBILITY.md) for the scope of each computation.

## Literature and novelty note

The starting conjectures are Alexander Berkovich and Aritram Dhar, *New Borwein-Type Conjectures*, arXiv:2407.13788, published in *Experimental Mathematics*.  Their paper explicitly proposes the five modulus-3 and two modulus-5 cases considered here.

The main analytic precedents are Chen Wang and Christian Krattenthaler's *An asymptotic approach to Borwein-type sign pattern theorems* (arXiv:2201.12415) and Liuquan Wang's *Sign Changes of Coefficients of Powers of the Infinite Borwein Product* (arXiv:2108.03932).

A targeted literature search dated 18 September 2026 did not find a prior proof of the six eventual theorems stated here.  That search is evidence, not a guarantee of priority.  See [`NOVELTY_NOTES.md`](NOVELTY_NOTES.md).

## AI-use disclosure

The mathematical exploration, manuscript drafting, code generation, execution, and internal checking were performed with OpenAI models.  The public audit in this repository is also AI-based.  The repository is intentionally explicit about this provenance so readers can weigh the result accordingly.

## References

See [`REFERENCES.md`](REFERENCES.md).
