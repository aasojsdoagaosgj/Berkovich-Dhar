# Infinite products for the six remaining cases

17 September 2026. Anonymous hobby research, derived and checked with
OpenAI Codex. No independent review or formal verification is claimed.

Write G_p(q)=prod_{l>=1,p not dividing l}(1-q^l). The following signs
hold at every degree, with zero admitted in either sign:

| p | power s | signs by residue 0,...,p-1 | zeros |
|---|---|---|---|
| 3 | 5,6,7,8 | +,-,+ | none |
| 5 | 2 | +,-,-,+,+ | degree 9 only |
| 5 | 3 | +,-,0,+,0 | exactly residues 2 and 4 |

Here is a proof using an explicit error theorem and a finite exact
certificate, rather than the sign table in the source.

## The published estimate and a common cutoff

Use [Liuquan Wang, arXiv:2108.03932v3, Theorem 1.1](https://arxiv.org/pdf/2108.03932).
All six cases satisfy s(p-1)<=24. Put mu=s(p-1)/24, d=j-mu and
X=4*pi*sqrt(mu*d)/p. For prime p, the theorem has A=1 and
M=sqrt(s(p-1))/(2p), and reads

```math
[q^j]G_p^s=\frac{2\pi\sqrt\mu}{p\sqrt d}\,C_{p,s}(j)I_1(X)+E,
\quad C_{p,s}(j)=\sum_{h=1}^{p-1}
 e^{-s\pi i\,s(h,p)-2\pi i jh/p}.
```

The Dedekind sum is denoted s(h,p) only in this formula; the exponent
s elsewhere is an integer power. The explicit error in that theorem is

```math
|E|\le \frac{\pi^{7/4}\mu^{1/4}}{2^{3/4}}e^{X/2}
 +2p e^{2+8\pi\mu}
 +2e^2p^{s/2+1}\exp\left\{\pi\mu+s\left(
 \frac{e^{-\pi}}{(1-e^{-\pi})^2}
 +\frac{e^{-\pi/p}}{(1-e^{-\pi/p})^2}\right)\right\}.
```

In our cases 1/3<=mu<=2/3, p<=5, s<=8. Using 31/10<pi<22/7,
e<3, e^{-pi}<1/20 and e^{-pi/5}<27/50 gives the common bound

```math
|E|<5e^{X/2}+10^{18}.                                      (I1)
```

Indeed the two constant terms are at most
10*3^19 + 2*9*5^5*3^25 < 10^18. The sum in parentheses is below
20/361+1350/529<8/3. All these rational comparisons are checked by
`certify_infinite.py`; bounds on exponentials use positive Taylor sums.

For p=3, s(1,3)=1/18 and s(2,3)=-1/18, so
C=2*cos(s*pi/18+2*pi*j/3). For p=5 the Dedekind sums are
(1/5,0,0,-1/5), so
C=2*cos(s*pi/5+2*pi*j/5)+2*cos(4*pi*j/5).
Elementary cosine values give the signs in the table. Every nonzero
C has |C|>1/3. In the modulus-3 cases the smallest absolute value is
2*sin(pi/18)>1/3; in the modulus-5 cases it is
(3-sqrt(5))/2>1/3. The only vanishing values are p=5,s=3,j=2,4 mod 5.

For X>=1 the Bessel integral gives I_1(X)>=e^X/(25*sqrt(X)):
restrict its positive integral to [1-1/X,1-1/(2X)], as in
`../quartic/INFINITE_SIGNS.md`. Therefore the magnitude of a nonzero
main term exceeds e^X/(300*d^{3/4}). The prefactor after substitution
is sqrt(pi)*mu^{1/4}/(25*sqrt(p)) times |C|; the stated bound follows
already from pi>3, mu>=1/3, p<=5 and |C|>1/3.

For d>=2000, X>(7/5)*sqrt(d), and the two upper bounds for the
error/main ratio are decreasing in d. At d=2000 they sum to less than

```math
450000e^{-31}+9\cdot10^{22}e^{-62}<1/100.                 (I2)
```

Here 2000^{3/4}<300 and (7/5)*sqrt(2000)>62. The inequality
e^{31}>2*10^{13}, certified by 100 Taylor terms, proves (I2).
Thus every j>=2001 with a nonzero C has the asserted strict sign.
The code verifies all coefficients through degree 2000 by exact
integer multiplication, and finds exactly the exceptional zero listed.

## The two identically zero residue classes

Jacobi's identity
(q;q)_infinity^3=sum_{k>=0}(-1)^k(2k+1)q^{k(k+1)/2}
has degrees congruent only to 0,1,3 modulo 5. Multiplication by
1/(q^5;q^5)_infinity^3 preserves residues. Therefore the coefficients
of G_5^3 in residues 2 and 4 are identically zero. This is an identity,
not an inference from finite zero samples or a vanishing main term.

The finite polynomial P_{n,p}^s agrees with G_p^s through degree pn,
since its first omitted factor is of degree pn+1. This completes the
initial-degree input needed for all six finite-polynomial problems.
