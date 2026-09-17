# Uniform growing-saddle estimates, including the vanishing cube residues

18 September 2026. Derived and written with OpenAI Codex. Ordinary proof
draft; no independent review. Notation is from MOD5_COMPACT.md.

The argument below applies to s=2,3. The weak case means s=3 and
a=2 or 4; the other cases are called nonvanishing. This distinction
is solely about the infinite product at the endpoint.

## 1. Exact tail and uniform expansion

Write G=G_5, P=P_{n,5}, T_n=P/G, w=z/(5n), x=e^{-z}, and

```math
C=2\pi^2/15,\qquad
R(z)=\frac{C-\operatorname{Li}_2(e^{-z})+
                      \operatorname{Li}_2(e^{-5z})/5}{z},\quad H=sR.
```

At a primitive root zeta^h the eta transformation gives

```math
G(\zeta^h e^{-w})=\kappa_h e^{A/w-w/6}U_h(w),\quad
A=2\pi^2/75,\quad
(\kappa_1,\kappa_2,\kappa_3,\kappa_4)=(e^{-i\pi/5},1,1,e^{i\pi/5}),
```

where U_h is the dual product in earlier manuscript (5.4). The
Dedekind sums are (1/5,0,0,-1/5). For w=w_0+iv, |v|<=3w_0/4,
that formula gives U_h=1+O(exp(-1/w_0)) uniformly as w_0 tends to
zero. Indeed the dual nome has modulus at most
exp(-(4*pi^2/25)*(16/25)/w_0), and its exponent exceeds 1/w_0.

The absolutely convergent logarithm of the inverse tail is exactly

```math
\log T_n(\zeta^h e^{-w})^s
=s\sum_{k\ge1}\frac{x^k}{k}
 \left\{\frac1{1-\zeta^{hk}e^{-kw}}-
                      \frac1{1-e^{-5kw}}\right\}.       (E1)
```

Set a_k=-1 if 5 does not divide k and a_k=4 otherwise. In the
closed sector |Im v|<=3 Re v/4, the bracket in (E1), with v=kw,
equals a_k/(5v)+b_{h,k}+O(|v|), uniformly for every v and k,
where b_{h,k}=1/(1-zeta^{hk})-1/2 for 5 not dividing k, and 0
otherwise. To see uniformity, divide the remainder by v: its
singularity at zero is removable by the Taylor series, its
denominators have no zeros in the sector away from zero, and for
large |v| the quotient is bounded. There are only five possible
root factors. The same statement holds in any fixed narrower
complex neighborhood needed below.

Consequently, for Re z=tau>=5 and |Im z|<=3tau/4,

```math
\log T_n(\zeta^h e^{-w})^s
=\frac{s}{5w}[-\operatorname{Li}_2(x)+\operatorname{Li}_2(x^5)/5]
 +s\lambda_h(z)-s\log\kappa_h+E_h(z),\quad
 E_h(z)=O(|w||x|).                                    (E2)
```

The phase identity in (E2) follows either by expansion of lambda_h,
or from sum_{b=1}^4 B_1(b/5)*zeta^{hbk}
=-1/(1-zeta^{hk})+1/2 when 5 does not divide k, and zero otherwise.
At x=0 the constant is log(kappa_h). The common radial exponent
after combining (E2) with the modular term is exactly nH(z).
Thus, up to the dual error, the root product is

```math
e^{nH(z)-s w/6}e^{s\lambda_h(z)+E_h(z)}.               (E3)
```

For the nonvanishing cases let B_{n,a}(z) be the sum over h of the
last two amplitude factors in (E3), weighted by zeta^{-ah}.
Then B_{n,a}=Psi_{s,a}+O(w) uniformly.

For a weak case define instead

```math
B_{n,a}(z)=e^z e^{-s w/6}
                \sum_{h=1}^4\zeta^{-ah}e^{s\lambda_h(z)+E_h(z)}.
```

Since sum_h zeta^{-ah}kappa_h^3=0 and E_h=O(w*x), this amplitude
is Psi_{3,a}/x+O(w), uniformly even when x tends to zero.
Its complex analyticity and boundedness hold on disks
|z-tau|<=tau/4. Those disks have Re z>=15/4; the same estimates
hold there with fixed enlarged constants. Cauchy's estimate gives
amplitude variation O(|z-tau|/tau) on smaller disks. Both weak
limiting amplitudes extend analytically to x=0.

On 0<=x<=e^{-5}, each nonvanishing Psi and each weak Psi/x has
the required strict sign and a nonzero minimum absolute value.
This follows from the factored formula (P1), the locations of its
simple zeros below tau=5, and the nonzero endpoint values or
derivatives. By compactness the minimum is a positive constant.
For sufficiently small w_0 the same uniform sign margins hold
for the real B_{n,a}.

## 2. Preserving the weak factor in every error

For a weak residue, [q^j]G^3=0 identically, so use

```math
[q^j]P^3=[q^j]\{G^3(T_n^3-1)\}                       (E4)
```

on the entire coefficient circle. The principal modular term of
the subtracted G^3 cancels exactly when the four root arcs are
summed with the same z, because sum_h zeta^{-ah}kappa_h^3=0.
The leading remaining root sum is therefore e^{nH(z)-z}B_{n,a}(z).
The dual errors, however, are bounded as
(U_h^3-1)*(T_n^3-1), before doing that sum. In particular they
retain the factor x_0=e^{-tau}.

On the whole circle |q|=e^{-w_0}, the positive tail majorant gives

```math
|\log T_n(q)^s|\le \frac{s}{w_0}\operatorname{Li}_2(x_0),
\quad |T_n(q)^s-1|\le
 \frac{s\operatorname{Li}_2(x_0)}{w_0}
             \exp\{s\operatorname{Li}_2(x_0)/w_0\}.     (E5)
```

Here sum_{l>=1,5 not dividing l} exp(-k*l*w_0)
<=1/(exp(k*w_0)-1)<=1/(k*w_0). Also
Li_2(x_0)<=x_0/(1-x_0)<=1/139 since e^5>140.
The common radial tail exponent in (E2) is at least
-s*Li_2(x_0)/(5*w_0). Therefore, after division by the radial
main term and, in the weak case, by x_0, the dual error is bounded
by C*w_0^{-1}*exp(-c/w_0), with c>0 independent of tau>=5.
For example (6s/5)*Li_2(x_0)<1/20 for s<=3, whereas the dual
exponent is at least 1. This justifies the weak-factor retention
even if tau grows as fast as a constant times sqrt(n).

## 3. The complementary circle

The dissection G(q)=g(q^5)^2-q*g(q^5)*h(q^5)-q^2*h(q^5)^2
and the positive-product angular estimate of earlier manuscript
Section 5.1 show that, outside all five arcs
|theta-2*pi*h/5|<=3w_0/4,

```math
|G(e^{-w_0+i\theta})|\le C_0
                  \exp\{A/w_0-1/(25w_0)\}.             (E6)
```

More explicitly the positive products are g^2,gh,h^2 at radius
exp(-5w_0). Each has loss at least 1/(5*(5w_0)); their real-axis
values are bounded by 2*exp(A/w_0) for small w_0 by earlier (5.7).
The three dissection terms give (E6). This is a uniform full-angle
estimate, not a sampled bound.

The arc at h=0, which is not a principal arc, is handled separately:
eta transformation gives G(e^{-w})=sqrt(5)*exp(-5A/w-w/6)U_0(w).
For |Im w|<=3w_0/4, Re(1/w)>=16/(25w_0), giving an exponentially
smaller modulus than the primitive-root main term.

Multiply (E6) by the tail majorant in (E5), using its second form
and (E4) in the weak cases. Relative to exp(nH(tau)), with x_0
also removed in the weak cases, the exponent loss is at least

```math
\frac{s}{w_0}\left(\frac1{25}-\frac65\operatorname{Li}_2(x_0)\right)
 >\frac{s}{40w_0},                                    (E7)
```

because 1/25-6/(5*139)>1/40. Prefactors are at worst C/w_0.
The h=0 arc has a larger loss by the preceding eta formula.
Thus the entire complementary circle retains the weak x_0 factor
and is exponentially negligible with parameter N=n/tau=1/(5w_0).

## 4. Gaussian control on the growing arcs

Write nH(z)=s*C*n/z+K(z), where
K=(sn/z)*sum_{k>=1}a_k*exp(-kz)/k^2. On Re z=tau>=5,

```math
|K''(z)|\le\frac{4sn}{\tau^3}
 (\tau^2+2\tau+2)\frac{e^{-\tau}}{1-e^{-\tau}}
 <\frac{11sn}{10\tau^3}.                              (E8)
```

For the inequality, the polynomial times exp(-tau) decreases for
tau>=5, and its value is bounded by 37/139. Similar termwise
differentiation gives fixed bounds for |(nH)'''|*tau^4/n and
|(nH)''''|*tau^5/n. With 13/10<C<4/3 this implies

```math
\frac{3sn}{2\tau^3}<g_n:=nH''(\tau)<\frac{4sn}{\tau^3}.
```

The exact reciprocal term and Taylor's integral remainder for K
give, for |y/tau|<=3/4,

```math
\Re\{nH(\tau-iy)-nH(\tau)\}
\le-sn(y/\tau)^2/\tau\,[16C/25-11/20]
\le-(s/4)N(y/\tau)^2.                                (E9)
```

The K' term is real on the real axis and contributes no real part.

Let d=j in a nonvanishing case and d=j-5n in a weak case. Choose
tau>=5 by -nH'(tau)=d/(5n). On the root arcs the extra factor
exp(-z) has already been absorbed into d in the weak case. Scaling
y=tau*v/sqrt(N), the quadratic exponent has fixed positive upper
and lower bounds; the cubic error is O(|v|^3/sqrt(N)); the bounded
analytic amplitudes vary by O(|v|/sqrt(N)). Gaussian integration,
(E9), and the complementary bounds prove uniformly

```math
\frac{c_j(n)}{E_{n,j}}=B_{n,a}(\tau)+O(N^{-1/2}),\quad
E_{n,j}=\frac{\exp\{nH(\tau)+d\tau/(5n)\}}
                         {5n\sqrt{2\pi g_n}}>0.        (E10)
```

For clarity, (E10) includes both the O(w_0) amplitude perturbation
and the dual errors described above. The latter are exponentially
small even after the O(N) weak prefactor. The complementary-circle
normalization contributes n*sqrt(g_n)=O(N^{3/2}), so no unbounded
power of n independent of N has been discarded.

## 5. Coverage

The saddle map -nH' is strictly decreasing to zero. For every
5n<=d<=-5n^2H'(5) it has a unique tau>=5. Moreover -H'(tau)
<s*C/tau^2. One way to see this last inequality is to write the
radial tail as -sn/tau times the positive decreasing integral
integral_tau^infinity log(1+e^{-t}+...+e^{-4t}) dt; its derivative
is positive. Thus the saddle equation implies tau^2<s*C*n<4n
and N>sqrt(n)/2. Hence (E10) fixes all signs uniformly for large n.

In nonvanishing cases this covers degrees from 5n to
-5n^2H'(5); initial degrees through 5n are covered by INFINITE_SIGNS.
In weak cases it covers degrees from 10n to 5n-5n^2H'(5);
MOD5_FIRST_TAIL covers all needed weak degrees below 10n, and
the degrees through 5n vanish exactly. Finally
-H'(5)>-H'(11/2), so both endpoint upper boundaries overlap the
compact range of MOD5_COMPACT by a positive multiple of n^2.
All estimates are uniform, so one sufficiently large n covers
every degree in the lower half with no integer gap.
