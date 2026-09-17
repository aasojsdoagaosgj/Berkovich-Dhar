# The two modulus-5 compact saddle theorems

18 September 2026. Anonymous research with OpenAI Codex. This document
uses the whole-circle localization proved in the existing manuscript
`../third-borwein-mod5/manuscript.md`, Section 8; that dependency must be retained
in any review package. This is not an independently reviewed result.

Let zeta=exp(2*pi*i/5), s=2 or 3, and

```math
R(z)=\int_0^1\log(1+e^{-zu}+e^{-2zu}+e^{-3zu}+e^{-4zu})\,du,
\quad H=sR,\quad\alpha=-H'.
```

The variance formula on {0,1,2,3,4} shows R''>0, R'(0)=-1 and
R'(infinity)=0. A lower-half degree j has unique saddle
alpha(tau)=j/(5n^2). The compact range here is 0<=tau<=11/2.

## Four roots and their combined sign

For h=1,...,4 and b=1,...,4, set
f_{h,b}(u,z)=log(1-zeta^{hb}exp(-zu)), continued from z=0. Put

```math
\lambda_h(z)=\sum_{b=1}^4 B_1(b/5)
                         (f_{h,b}(1,z)-f_{h,b}(0,z)),
\quad \Psi_{s,a}(z)=\sum_{h=1}^4\zeta^{-ah}e^{s\lambda_h(z)}.
```

On the real axis write x=e^{-tau},
A=atan(sin(2*pi/5)*x/(1-cos(2*pi/5)*x)),
B=atan(sin(4*pi/5)*x/(1-cos(4*pi/5)*x)), and

```math
U=-\pi/10+(A+2B)/5,\quad V=-\pi/10+(2A-B)/5.
```

Pairing conjugate roots and factoring the sum of two cosines gives

```math
\boxed{\Psi_{s,a}=4\cos(sU-3\pi a/5)\cos(sV+\pi a/5).} (P1)
```

Both U,V increase strictly with x from -pi/10 at x=0 to 0 at x=1.
For V, use A'=sin(2*pi/5)/(1-2*cos(2*pi/5)*x+x^2), the analogous
B', and 2*A'>B'. The latter follows from the larger numerator and
smaller positive denominator of A'. Thus U,V decrease strictly with
tau. Their endpoint derivatives with respect to x are nonzero.

For s=2,3, (P1) immediately gives Psi_0>0, Psi_1<0, Psi_2<0
for 0<x<=1. The single zeros in the other two residues are

```math
V(e^{-\tau_{s,3}})=-\pi/(10s),\qquad
U(e^{-\tau_{s,4}})=-\pi/(10s).                          (P2)
```

For residue 3 the first cosine is positive and the second changes
once; for residue 4 the second cosine is negative and the first
changes once. Each Psi has positive tau-derivative at its zero.
Both zeros lie in (0,5), and both amplitudes are positive for tau
above their zero and negative below it. One can see the bound 5
without a decimal root calculation: at x<=1/140, A,B<1/100, so
U,V are within 1/100 of -pi/10 and are smaller than -pi/20;
at tau=0 they are zero. Here e^5>140 follows from its Taylor sum.

At x=0 the s=2 amplitudes have the signs +,-,-,+,+ with nonzero
margins. At s=3 the amplitudes in residues 2 and 4 vanish, but
Psi_{3,2}/x and Psi_{3,4}/x have finite nonzero limits with signs
negative and positive respectively. In (P1) precisely one cosine
factor vanishes at that endpoint, with nonzero x derivative.

## Uniform expansion to order n^{-3}

Shifted Euler--Maclaurin at each root gives

```math
\log P_{n,5}(\zeta^h e^{-z/(5n)})^s
 =nH+s\lambda_h+b_h/n+d_h/n^2+O(n^{-3}),              (P3)
```

where

```math
b_h=\frac s2\sum_{b=1}^4B_2(b/5)
 (\partial_u f_{h,b}(1,z)-\partial_u f_{h,b}(0,z)),
\quad d_h=\frac s6\sum_{b=1}^4B_3(b/5)
 (\partial_u^2 f_{h,b}(1,z)-\partial_u^2 f_{h,b}(0,z)).
```

All functions are analytic and their required derivatives uniformly
bounded in a small complex neighborhood of [0,11/2]. The shifted
Bernoulli-integral remainder proves the asserted uniformity.

Section 8 of the earlier modulus-5 manuscript proves, for n>=31147,

```math
\log|P_{n,5}(e^{-\tau/(5n)+i\theta})|
 \le n(R(\tau)-1/1000)                              (P4)
```

outside the four arcs |theta-2*pi*h/5|<= (6/5)/(5n), uniformly
on this tau interval. This is a bound on the base product, so it
may be raised to either power. Its proof uses a full Dirichlet
cover and Fourier smoothing, with interval certificates; it is
not a finite grid of angle samples.

These arcs may be reduced to a small fixed width eta/(5n).
On their remaining closed portions eta<=|Im z|<=6/5, the strict
triangle inequality for 1+exp(-zu)+...+exp(-4zu), integrated over
0<u<1, gives Re R(tau-iy)<R(tau). Compactness gives a positive
uniform gap. All the Euler--Maclaurin factors stay nonzero here,
since 6/5<2*pi/5. Thus their O(1) amplitudes do not change the gap.

Put g=H'' and B_{n,j}=exp(n(H+alpha*tau))/(5n*sqrt(2*pi*n*g)).
The local Gaussian argument from quartic Section 3, now summing
four roots, yields

```math
c_j/B_{n,j}=\Psi_{s,a}+A_{1,s,a}/n+A_{2,s,a}/n^2+O(n^{-3}). (P5)
```

For an explicit specification put
B_a=sum_h zeta^{-ah}exp(s*lambda_h)*b_h and
D_a=sum_h zeta^{-ah}exp(s*lambda_h)*(d_h+b_h^2/2). Then
A_1=B_a+S_1(Psi), A_2=D_a+S_1(B_a)+S_2(Psi), where S_1,S_2
are the Gaussian differential operators in quartic (C7),(C7a),
with H=sR. In particular A_1,A_2 are real analytic with bounded
first derivatives. No phase-dependent b_h has been incorrectly
replaced by a common scalar.

The proof of the remainder uses the expansion in n^{-1/2} through
degree 5, cancellation of odd moments, and an integrable uniform
O(n^{-3}) error on the central Gaussian region. The tails of that
region and (P4) are exponentially small after normalization.

Consecutive degrees in one residue differ by 5, giving
Delta tau=-1/(n^2*H'')+O(n^{-4}). Subtracting (P5) at these two
degrees gives

```math
\Delta(c_j/B_{n,j})=-\Psi_{s,a}'/(n^2H'')+O(n^{-3}).     (P6)
```

Near each zero (P2) this is strictly negative for large n. Away
from these neighborhoods the signs are fixed by (P1). Hence in
the compact range each of residues 3,4 has exactly one transition,
allowing at most one zero at that transition, while residues 0,1,2
keep their signs. The shift index satisfies

```math
k_{n,s,a}=\alpha_{s,a}n^2+\beta_{s,a}n+O(1),\quad
\alpha_{s,a}=-H'(\tau_{s,a}),\quad
\beta_{s,a}=H''(\tau_{s,a})A_{1,s,a}(\tau_{s,a})/
                                      \Psi_{s,a}'(\tau_{s,a}). (P7)
```

The passage j=5k+a affects only the bounded term. Formula (P7)
does not by itself prove the signs outside this compact range;
the endpoint and initial-tail arguments supply that connection.
