# Uniform control of the small-degree saddle range

17 September 2026. Anonymous derivation with OpenAI Codex. This is an
ordinary analytic proof, not an independently reviewed or formal proof.

The notation $P_n,F_n,c_j$ is as in COMPACT_PROOF.md. In this document
$a_0=\pi^2/3$, and $a\in\{0,1,2\}$ always denotes a residue.

## 1. Exact modular and tail formulas

Let $w=z/(3n)$, $\Re z>0$, and $t(z)=e^{-4a_0n/z}$.
Dedekind eta transformation gives

```math
F_\infty(\omega e^{-w})
=e^{4\pi^2/(27w)-w/3-2\pi i/9}F_\infty(\omega^2e^{-4\pi^2/(9w)}).  (E1)
```

Here $F_\infty=P_\infty^4$. This follows directly, without taking
an unspecified cube root: put $\sigma=1/3+iw/(2\pi)$,
$\sigma'=\sigma/(1-3\sigma)=-1/3+2\pi i/(9w)$.
The eta ratio transforms by $e^{-i\pi/6}$. Since
$P_\infty(e^{2\pi i\sigma})=e^{\pi i\sigma/6}\eta(\sigma)/\eta(3\sigma)$,
the resulting multiplier for $P_\infty$ is
$e^{\pi^2/(27w)-w/12-i\pi/18}$. Raise this identity to the fourth
power. The eta multiplier uses $s(-1,3)=-1/18$ for the matrix
$(1,0;3,1)$, with the principal square root in
[DLMF 23.18.5--7](https://dlmf.nist.gov/23.18#E5).
The other root follows first on the real axis by conjugation, then
as a separate analytic identity by continuation.

Define $a_k=1$ for $3\mid k$ and $a_k=-1/2$ otherwise,
and $\chi_3(k)=1,-1,0$ for residues $1,2,0$. Put

```math
K_n(z)=4\sum_{k\ge1}\frac{a_k e^{-kz}}k
          \sum_{l\ge1,\,3\nmid l}e^{-klz/(3n)},
\quad V_n(z)=2\sqrt3\sum_{k\ge1}
 \frac{\chi_3(k)e^{-kz}}{k(1+2\cosh(kz/(3n)))},
```

```math
L_n(z)=\frac{4a_0n}{3z}-\frac z{9n}+K_n(z),\qquad
T_a(z)=2\cos(-2\pi/9-2\pi a/3+V_n(z)).
```

The logarithms of the omitted-factor inverse tails at the two roots are
exactly $K_n\pm iV_n$, by their absolutely convergent series. Combining
them with (E1), the weighted sum of the two finite root products is

```math
\omega^{-a}F_n(\omega e^{-z/(3n)})+
\omega^aF_n(\omega^2e^{-z/(3n)})
=e^{L_n(z)}\bigl(T_a(z)+O(e^{-cN})\bigr)               (E2)
```

uniformly for $z=\tau-iy$, $\tau\ge6$, $|y|\le3\tau/2$,
and $N=n/\tau\to\infty$. Indeed,
$|t(z)|\le e^{-16a_0N/13}$; the dual infinite products differ from
1 by $O(|t|)$ uniformly, by their convergent products. The phase
factors in (E2) are uniformly bounded, as shown below. All constants
in this document are absolute and independent of $\tau\ge6$.

## 2. Exponent and amplitude bounds

Differentiating the absolutely convergent series and comparing the
inner sum to an integral gives, for $0\le p\le4$,

```math
|K_n^{(p)}(\tau-iy)|\le12n\sum_{l=0}^p
 \frac{p!}{(p-l)!\tau^{l+1}}
 \sum_{k\ge1} k^{p-l-2}e^{-k\tau}.                    (E3)
```

For completeness, enlarge the sum over positive $l$ to all positive
integers, and compare $(1+s)^p e^{-k\tau s}$ to its integral on
$[0,\infty)$. It is decreasing because $k\tau\ge6>p$.
Expanding the polynomial before integrating gives (E3).
On the real axis one can also write

```math
K_n(\tau)=-2\sum_{l>3n,\,3\nmid l}
 \log(1+e^{-l\tau/(3n)}+e^{-2l\tau/(3n)}).
```

Thus $K_n<0,K_n'>0,K_n''<0$. Write $g=L_n''(\tau)$.
Since $x=e^{-\tau}<1/400$, (E3) implies

```math
|K_n''|\le\frac{600}{399}\frac n{\tau^3},\qquad
\frac{20}{3}\frac n{\tau^3}<g<\frac{28}{3}\frac n{\tau^3},
\qquad |L_n'''|<39n/\tau^4.                          (E4)
```

These elementary margins use $39/4<\pi^2<10$. For the third
derivative, use $\sum kx^k=x/(1-x)^2$, reduce each
$\tau^j e^{-\tau}$ to its value at 6, and bound the reciprocal
term separately. The same argument gives a fixed bound for
$|L_n''''|\tau^5/n$. In particular every scaled quantity needed
for a third-order local Taylor estimate is uniformly bounded.

For $u=y/\tau$, $|u|\le3/2$, the exact reciprocal term and
the second-derivative bound for $K_n$ give

```math
\Re(L_n(\tau-iy)-L_n(\tau))
\le-Nu^2\left(\frac{4a_0}{3(1+u^2)}-\frac{300}{399}\right)
\le-\frac8{15}Nu^2.                                 (E5)
```

The linear term $-z/(9n)$ has no real-part change. The first
derivative of $K_n$ on the real axis is real, so its linear term
also contributes no real part.

For $s=A+iD$ with $A>0, |D|\le3A/2$,
$|1+2\cosh s|\ge1$. If $A\le1$, use its real part and
$\cos(3/2)>0$; if $A\ge1$, use
$2|\cosh s|-1\ge2\sinh A-1>1$.
It follows that $V_n$, and hence $T_a$, are uniformly bounded
both on the root arcs and on $|z-\tau|\le\tau/2$, where
$\Re z\ge3$. Cauchy's estimate therefore gives
$|T_a'(\tau)|\le C/\tau$, and the same bound for its variation
on a smaller disk.

On the real axis the nonzero terms of the phase series alternate with
decreasing magnitude. Consequently

```math
0<V_n(\tau)\le\frac{2\sqrt3 x}{1+2\cosh(\tau/(3n))}
 \le\frac{2\sqrt3}{3}x<1/100.
```

This yields uniform sign margins

```math
T_0(\tau)>1,\qquad T_1(\tau)<-1,\qquad T_2(\tau)>1/4.       (E6)
```

For the last inequality use
$T_2=2\sin(\pi/18-V_n)$, $\pi>3$, and
$\sin s\ge s-s^3/6$. No small cancelling amplitude occurs here.

## 3. The rest of the coefficient circle

Use (C3)--(C4) of COMPACT_PROOF.md with
$r=e^{-1/(3N)}$. Take the two root arcs of half-width
$1/(2N)=\tau/(2n)$, exactly the arcs in (E2).
On the complement where $|1-re^{i\theta}|<1/3$, (C4) gives
$D\ge cN-C$, since $1-r^{3n}=1-e^{-\tau}\ge1-e^{-6}$.
Elsewhere the negative term in (C3) is bounded by 2.4.
If the nearest root is nontrivial, then

```math
|\rho|\ge\frac{1-e^{-6}}{2N(1-e^{-1/N})}>\frac{1-e^{-6}}2.
```

If the nearest root is 1, the condition $|1-re^{i\theta}|\ge1/3$
keeps its angular distance away from zero when $N$ is large;
the same lower bound follows. Finally
$(1-r^{n/2})/(1-r^3)\ge(1-e^{-1})N$.
Thus again $D\ge cN-C$. Raising to the fourth power shows that
the modulus of the full finite product outside the two arcs is at most
$C e^{-cN}|F_n(\omega r)|$.
By (E1), $|F_n(\omega r)|/e^{L_n(\tau)}=1+O(e^{-cN})$.
The entire complementary integral is therefore exponentially small
relative to the radial main exponential.

## 4. Uniform endpoint coefficient estimate

Suppose $\tau\ge6$ solves
$-L_n'(\tau)=j/(3n)$, and $N=n/\tau\to\infty$.
For $j\equiv a\pmod3$, Cauchy's coefficient integral, (E2), and
the saddle equation give

```math
\frac{c_j(n)}{\mathcal E_{n,j}}
=T_a(\tau)+O(N^{-1/2}),\qquad
\mathcal E_{n,j}=
 \frac{e^{L_n(\tau)+j\tau/(3n)}}{3n\sqrt{2\pi g}}>0.       (E7)
```

Here the error is uniform in the whole stated range, including growing
$\tau$. To see this explicitly, scale $y=\tau v/\sqrt N$.
The quadratic exponent has a coefficient bounded above and below by
positive constants by (E4); the cubic remainder is
$O(|v|^3/\sqrt N)$; and the amplitude variation is
$O(|v|/\sqrt N)$. On $|v|\le N^{1/20}$ their Gaussian-weighted
integrals are $O(N^{-1/2})$. On the rest of the arcs (E5) gives
a Gaussian tail. Dual errors contribute $O(\sqrt N e^{-cN})$.
The complementary circle from Section 3 contributes
$O(n\sqrt g\,e^{-cN})=O(N^{3/2}e^{-cN})$, using (E4).
All are uniform and tend to zero. This proves (E7), so (E6) determines
the coefficient signs whenever $N$ is sufficiently large.

## 5. Coverage and overlap with the compact range

The function $-L_n'$ is strictly decreasing on $[6,\infty)$,
with limit $1/(9n)$. Consequently the endpoint saddle exists uniquely
whenever

```math
3n\le j\le j_e(n):=-3nL_n'(6).
```

Because $K_n'>0$, its equation implies

```math
1\le\frac j{3n}
<\frac{4a_0n}{3\tau^2}+\frac1{9n},\qquad
\tau^2<5n,\qquad N>\sqrt{n/5}.                       (E8)
```

Thus (E7) gives the required $+,-,+$ signs on this entire interval
for all sufficiently large $n$.

At the fixed argument 6, Euler--Maclaurin from COMPACT_PROOF.md, or its
real part and (E1), gives
$L_n'(6)=nH'(6)+O(n^{-1})$. The derivative of the dual logarithm
is exponentially small. It follows that

```math
\frac{j_e(n)}{3n^2}=\alpha(6)+O(n^{-2})>\alpha(7)
```

for all sufficiently large $n$, because $H''>0$.
The compact range begins at $3n^2\alpha(7)$. The overlap has a
positive limiting width on the $n^2$ scale, so no integer rounding
gap arises. Together the two ranges cover every degree from $3n$
to the midpoint $6n^2$. Initial degrees are covered by INFINITE_SIGNS.md.
