# Four eventual sign-pattern theorems, modulus 3 and powers 5 through 8

18 September 2026. Anonymous hobby research. Derivation, manuscript,
code and checks by OpenAI Codex; independent review and formalization
have not been performed. This is a conventional proof draft.

## Statement

For s in {5,6,7,8}, let

```math
F_{n,s}(q)=\prod_{l=1,\,3\nmid l}^{3n}(1-q^l)^s
 =\sum_{j=0}^{3sn^2}c_j^{(s)}(n)q^j.
```

There exists an integer N_s such that for every n>=N_s:

* all coefficients c_{3k}^{(s)}(n) are strictly positive;
* in residue 2 there is a single transition from positive to negative,
  allowing at most one zero at the transition;
* all coefficients in residue 1 are determined by reciprocity.

Let k_{n,s} be the last nonnegative residue-2 index. Then

```math
k_{n,s}=\alpha_s n^2+\beta_s n+O(1),\quad
\alpha_s=-sR'(\tau_s),\quad
\tau_s=-\log\frac{2\tan(\pi/6-\pi/(2s))}
                         {\sqrt3-\tan(\pi/6-\pi/(2s))},
\quad R(z)=\int_0^1\log(1+e^{-zu}+e^{-2zu})\,du.       (M1)
```

This establishes the qualitative and convergence assertions of
[Berkovich--Dhar, Conjecture 2.1](https://arxiv.org/html/2407.13788v3#S2)
for these powers and sufficiently large n. It does not assert all n,
and no numerical value of N_s is supplied. Numerical evaluations of
the exact constants (M1), not interval certificates, are

| s | alpha_s |
|---|---|
| 5 | 1.2692836849298353885... |
| 6 | 1.7772890339199816431... |
| 7 | 2.2814202796053000080... |
| 8 | 2.7838683312383998477... |

The first three printed decimals in the conjecture differ from these
limits. The proof concerns the exact expressions, not the printed
decimal estimates.

## Compact saddles, including actual adjacent coefficient comparisons

We use the analytic expansion proved in `../quartic/COMPACT_PROOF.md`,
with the following substitutions, each made at the logarithm level:

```math
H=sR,\quad h_s=(s/4)h_4,\quad b_s=(s/4)b_4,\quad d_s=(s/4)d_4.
```

The Euler--Maclaurin formula for log(F_{n,s}) is therefore
nH +/- h_s + b_s/n +/- d_s/n^2 + O(n^{-3}), uniformly in a
complex neighborhood of [0,7]. The same derivative bounds follow
from the same fixed compact set. Set

```math
\phi_s(\tau)=-s\pi/18+(s/3)\arctan
 \frac{\sqrt3e^{-\tau}}{2+e^{-\tau}},\qquad
\Psi_{s,a}=2\cos(\phi_s-2\pi a/3).
```

The phase decreases from 0 to -s*pi/18 and
phi_s'=-s*sqrt(3)*x/(6*(1+x+x^2))<0. Hence Psi_{s,0}>0
and Psi_{s,1}<0 on the whole nonnegative real axis. Psi_{s,2}
has exactly one zero, at phi_s=-pi/6, with positive derivative;
solving that equation gives (M1). This zero lies between 0 and 6:
the target tangent is between tan(pi/15) and tan(5*pi/48), so
e^{-tau_s} is bounded away from both 0 and 1.

The full-circle modulus estimates in quartic equations (C3),(C4)
are [Wang--Krattenthaler's Lemmas 9.2 and 9.4](https://arxiv.org/html/2201.12415#S9)
for the base product itself. Raising their modulus bounds to power s
therefore proves the same exponential localization. No assertion
restricted to powers 1,2,3 is used here.

Put alpha=-H', g=H''>0. For the saddle alpha(tau)=j/(3n^2),
the Gaussian normalization is

```math
B_{n,j}=\frac{e^{n(H(\tau)+\alpha\tau)}}
                   {3n\sqrt{2\pi n g(\tau)}}>0.
```

The local integration argument in quartic Section 3 now gives

```math
c_j/B_{n,j}=\Psi_{s,a}+A_{1,s,a}/n+A_{2,s,a}/n^2+O(n^{-3}), (M2)
```

uniformly for tau in [0,7]. A_1 and A_2 are the functions in
quartic (C7),(C7a) with precisely the substitutions listed above.
They are real analytic, with bounded first derivatives on this
interval. The proof of (M2) uses Taylor expansion through degree 5
in n^{-1/2}, vanishing odd Gaussian moments, a uniform integrable
O(n^{-3}) remainder, and the exponentially small complementary arcs.
Multiplication by the fixed s/4 does not alter any remainder order.

At two consecutive indices in a fixed residue class,

```math
\tau(j+3)-\tau(j)=-\frac1{n^2H''(\tau(j))}+O(n^{-4}).
```

Subtracting the two value expansions (M2), without differentiating
their remainder, gives

```math
\frac{c_{j+3}}{B_{n,j+3}}-\frac{c_j}{B_{n,j}}
=-\frac{\Psi_{s,a}'(\tau(j))}{n^2H''(\tau(j))}+O(n^{-3}). (M3)
```

In a fixed neighborhood of tau_s, Psi_{s,2}' is bounded below by
a positive constant. Thus (M3) is strictly negative for sufficiently
large n. Away from that neighborhood, the nonzero main-term margin
in (M2) fixes the sign. This proves exactly one transition throughout
the compact range and at most one zero there. At tau_s, solving
(M2) to order 1/n yields

```math
\beta_s=-\frac{\Psi_{s,2}''(\tau_s)}{2\Psi_{s,2}'(\tau_s)}
                  +\frac{H'''(\tau_s)}{2H''(\tau_s)}.       (M4)
```

Indeed A_1 at the zero consists only of its Psi' and Psi'' terms;
the terms with Psi vanish. The relation j=3k+2 changes only O(1).

## Growing saddles and initial degrees

All infinite coefficients have the strict signs +,-,+ by
`INFINITE_SIGNS.md`, including its finite exact certificate. Thus
the same signs hold for finite coefficients through degree 3n.

For tau>=6 use the exact modular-tail formulas of
`../quartic/ENDPOINT_PROOF.md`, with

```math
L_{n,s}=(s/4)L_{n,4},\quad V_{n,s}=(s/4)V_{n,4},\quad
T_{s,a}=2\cos(-s\pi/18-2\pi a/3+V_{n,s}).
```

These identities follow by multiplying the base logarithms by s;
in particular the eta multiplier is exp(-s*pi*i/18), and the exact
radial modular exponent is s*pi^2/(27*w)-s*w/12. The discarded dual
products still have a uniform O(exp(-c*N)) error, where N=n/tau.

All derivative and arc-decay bounds of quartic (E3)--(E5) scale by
s/4. Since s belongs to a fixed finite set, their constants remain
uniform. In particular L'' is bounded above and below by positive
constant multiples of n/tau^3 and the local real-part decay is
at least (s/4)*(8/15)*N*(y/tau)^2. The outside-arc argument applies
to the base product and hence to its s-th power. The prefactor after
normalization is O(N^{3/2}), not a separate growing power of n.

On the real axis 0<V_{n,s}<s/400<=1/50. The angle intervals then
give the uniform margins

```math
T_{s,0}>1/3,\qquad T_{s,1}<-1,\qquad T_{s,2}>1/2.       (M5)
```

For the first, the worst angle magnitude is 4*pi/9, and
2*cos(4*pi/9)=2*sin(pi/18)>1/3. For the second, the angles lie
between -10*pi/9 and -17*pi/18+1/50, inside (-4*pi/3,-2*pi/3).
For the third, angles after adding 2*pi lie between 2*pi/9 and
7*pi/18+1/50, inside (-pi/3,5*pi/12); cosine at 5*pi/12 exceeds 1/4.

For the saddle -L_{n,s}'(tau)=j/(3n), the same integral estimate gives

```math
\frac{c_j}{E_{n,j}}=T_{s,a}(\tau)+O(N^{-1/2}),\qquad
E_{n,j}=\frac{e^{L_{n,s}(\tau)+j\tau/(3n)}}
                   {3n\sqrt{2\pi L_{n,s}''(\tau)}}>0.      (M6)
```

The error is uniform for tau>=6. If 3n<=j<=-3nL_{n,s}'(6),
the saddle exists uniquely for all sufficiently large n, since
-L' decreases to s/(36n)<1. As K_{n,s}'>0, its equation gives
1<=s*pi^2*n/(9*tau^2)+s/(36n), whence tau^2<10n for large n.
Thus N>sqrt(n/10) uniformly, and (M5),(M6) fix all signs on
this range.

At the upper endpoint, Euler--Maclaurin gives
-L_{n,s}'(6)/n=-H'(6)+O(n^{-2}). Strict positivity of H''
implies -H'(6)>-H'(7), so this range overlaps the compact range
3n^2*(-H'(7))<=j<=3sn^2/2. No rounding gap is possible for large n.

## Completion for each of the four powers

The initial, growing-saddle and compact ranges cover the entire
lower half. There residue 0 is positive everywhere, residue 1
negative everywhere, and residue 2 changes once from positive to
negative. Since the total degree is 3sn^2 and the number of factors
is even, c_j=c_{3sn^2-j}. The upper-half residue 2 is therefore the
reflection of lower-half residue 1 and is negative; residue 0 stays
positive by reflection. This establishes the stated entire-polynomial
theorem. There are finitely many uniform sufficiently-large conditions
and four powers, so one common integer threshold also exists.

Numerical profile evaluations and finite-polynomial samples are
diagnostics only. The claimed proof is the chain of uniform analytic
estimates above together with the stated infinite-product certificate.
