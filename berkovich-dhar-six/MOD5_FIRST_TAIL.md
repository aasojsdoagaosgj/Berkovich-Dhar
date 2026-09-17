# The first nonzero tail for the modulus-5 cube

18 September 2026. Derived with OpenAI Codex; no independent review.
All inequalities between series below mean coefficientwise inequalities.

Let Q=q^5 and use the Rogers--Ramanujan series

```math
R_t(Q)=\sum_{k\ge0}\frac{Q^{k^2+(t-1)k}}{(Q;Q)_k},
\quad g=R_1,\ h=R_2,\ u=g-h=QR_3.
```

The recurrence R_t=R_{t+1}+Q^tR_{t+2} implies R_t>=R_{t+1}>=0.
The classical dissection, derived in Wang's Proposition 3.6 and in
`../third-borwein-mod5/manuscript.md` Section 2, is
G_5(q)=g(Q)^2-q*g(Q)*h(Q)-q^2*h(Q)^2.
Cubing this polynomial identity gives

```math
G_5(q)^3=A(Q)+qB(Q)+q^3C(Q),
\quad A=g^6-3Qgh^5,\ B=-3g^5h-Qh^6,\ C=5g^3h^3.       (F1)
```

In particular the omitted two residues vanish exactly.

## A bounded positive kernel

Put D=(gh+h^2-g^2)/(1-Q). The recurrence gives the positive identity

```math
D=\sum_{j\ge0}Q^{2j(j+1)}R_{2j+2}R_{2j+3}\ge0.         (F2)
```

For completeness this follows by setting
W_t=R_{t+1}^2-QR_tR_{t+2}; direct substitution gives
W_t=(1-Q)R_{t+1}R_{t+2}+Q^{2t+2}W_{t+2}.
Iterate from t=1. Its remainder tends to zero formally, since its
minimal degree tends to infinity; W_1=gh+h^2-g^2.

A second bound, needed for the cube, is D<=h^2. In fact

```math
h^2-D=\frac{Q(gR_3-h^2)}{1-Q}\ge0.                      (F3)
```

To prove the last inequality directly, symmetrize the two products:

```math
\frac{R_1R_3-R_2^2}{1-Q}
=\sum_{0\le j<k}\frac{Q^{j^2+k^2+2j}}{(Q;Q)_j(Q^2;Q)_{k-1}}
              (1+Q+\cdots+Q^{k-j-1})^2\ge0.
```

The removed factor (1-Q) from (Q;Q)_k combines with the outside
denominator to cancel the square (1-Q^{k-j})^2. This identity
establishes (F3) coefficientwise, not just on the positive real axis.

## The two partial-sum kernels

Define J_2=(A+B+QC)/(1-Q), J_4=(A+B+C)/(1-Q). Elementary
polynomial expansion, using g^2-gh-h^2=-(1-Q)D, gives

```math
J_4=3gh^5+h^6-(1-Q)^2D^3,                              (F4)
```

```math
-J_2=(1-Q)^2D^3+P,\quad
P=5g^3h^3-3gh^5-h^6
 =h^6+12uh^5+15u^2h^4+5u^3h^3.                        (F5)
```

Now J_4>=3gh^5-Q^2h^6>=2gh^5. The first inequality follows from
0<=D<=h^2 and -(1-Q)^2D^3>=-(1+Q^2)h^6. For the second,
Q^2h<=g: h<=g, and g has the product 1/(Q,Q^4;Q^5)_infinity,
whose coefficients are nondecreasing because part 1 is allowed.
Thus J_4 is strictly positive in every degree.

Also h=R_3+Q^2R_4 and h has part 2 available, so
h^6=R_3h^5+Q^2R_4h^5<=2R_3h^5. Consequently
2Qh^6<=4uh^5. Equation (F5) implies

```math
-J_2\ge P-2Qh^6
 \ge h^6+8uh^5+15u^2h^4+5u^3h^3>0.                    (F6)
```

The strict sign holds in every degree: h^6 has positive constant
coefficient and positive coefficients in every degree >=2 (parts 2
and 3 are available); uh^5=QR_3h^5 has positive degree-1 coefficient.

## Exact connection to the finite cube

Let T_n=prod_{l>5n,5 not dividing l}(1-q^l)^{-1}. Formally,

```math
T_n^3=1+3q^{5n}\frac{q+q^2+q^3+q^4}{1-q^5}+O(q^{10n+2}).
```

Since P_{n,5}^3=G_5^3*T_n^3, equations (F1) and this expansion give

```math
[q^{5n+5K+a}]P_{n,5}^3=3[Q^K]J_a,
\qquad a=2,4,\quad 0\le K\le n-1.                    (F7)
```

Thus residue 2 is strictly negative and residue 4 strictly positive
through this whole first tail. Before it, both are identically zero
through degree 5n. The formula includes the first nonzero indices
5n+2 and 5n+4 and extends through 10n-3 and 10n-1 respectively.
It is not asserted beyond its stated range.
