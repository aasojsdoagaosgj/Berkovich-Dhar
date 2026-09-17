# A uniform compact expansion and one sign change

17 September 2026. Anonymous research derivation written with OpenAI Codex.
The proof is ordinary mathematics; numerical samples are not its justification.

Write $P_n(q)=\prod_{j=1}^n(1-q^{3j-2})(1-q^{3j-1})$,
$F_n=P_n^4$, $c_j(n)=[q^j]F_n$, and $\omega=e^{2\pi i/3}$.
Set

```math
R(z)=\int_0^1\log(1+e^{-zu}+e^{-2zu})\,du,\quad H=4R,
\quad \alpha(\tau)=-H'(\tau).
```

The branches below are continued from real $z\ge0$ in a sufficiently
small complex neighborhood of $[0,7]$. All factors involved are
nonzero there. At zero the functions have their continuous values.

## 1. The real saddle and the two amplitudes

For $s\ge0$, let $J_s$ have probabilities proportional to
$1,e^{-s},e^{-2s}$ on $\{0,1,2\}$. Then

```math
R''(\tau)=\int_0^1u^2\operatorname{Var}(J_{\tau u})\,du>0,
\quad R'(0)=-1/2,\quad R'(\infty)=0.
```

Thus $\alpha$ decreases smoothly from 2 to 0. For every
$0<j/(3n^2)\le2$, there is a unique $\tau\ge0$ with
$\alpha(\tau)=j/(3n^2)$.

Let $f_b(u,z)=\log(1-\omega^b e^{-zu})$, for $b=1,2$, and

```math
h(z)=\frac23\{f_2(1,z)-f_2(0,z)-f_1(1,z)+f_1(0,z)\}.
```

The shifted Euler--Maclaurin formula for the sample positions
$(k+b/3)/n$, multiplied by 4, gives

```math
\log F_n(\omega^{\pm1}e^{-z/(3n)})
=nH(z)\ \pm h(z)+\frac{b(z)}n
 \ \pm\frac{d(z)}{n^2}+O(n^{-3}),                  (C1)
```

uniformly on a complex neighborhood of $[0,7]$, with uniform
derivative bounds of any fixed order on smaller neighborhoods. Here

```math
b(z)=-\frac z9\left(1-\frac{x(1+2x)}{1+x+x^2}\right),\qquad x=e^{-z}.
```

One may take
$d(z)=-2i\sqrt3 z^2x(1-x^2)/(81(1+x+x^2)^2)$.
To justify the asserted
uniformity, apply Euler--Maclaurin with a fixed number of terms to
$f_b(u,z)$. These functions and all the needed derivatives are bounded
on a compact complex neighborhood of $[0,1]\times[0,7]$. The usual
Bernoulli-integral remainder is therefore uniformly bounded. Further
terms can be retained before applying Cauchy's estimate for derivatives.
The common $1/n$ term uses $B_2(1/3)=B_2(2/3)=-1/18$;
the opposite $1/n^2$ terms use $B_3(1/3)=-B_3(2/3)$.

On the real axis $h=i\phi$, where

```math
\phi(\tau)=-\frac{2\pi}{9}+
 \frac43\arctan\frac{\sqrt3 e^{-\tau}}{2+e^{-\tau}},
\qquad \phi'(\tau)=-\frac{2\sqrt3 e^{-\tau}}{3(1+e^{-\tau}+e^{-2\tau})}<0.
```

For residue $a\in\{0,1,2\}$, define the analytic combined amplitude

```math
\Psi_a(z)=\omega^{-a}e^{h(z)}+\omega^a e^{-h(z)}.
```

For real $\tau$, this is $2\cos(\phi(\tau)-2\pi a/3)$.
In particular, $\Psi_0>0$, $\Psi_1<0$, and $\Psi_2$ has one
simple zero at $\tau_*$, with $\Psi_2'(\tau_*)>0$. Solving it gives

```math
x_*:=e^{-\tau_*}=\frac{2\tan(\pi/24)}{\sqrt3-\tan(\pi/24)},
\qquad \alpha_*=-4R'(\tau_*).                         (C2)
```

The zero lies in $(0,6)$. Below $\tau_*$, $\Psi_2<0$; above it,
$\Psi_2>0$.

## 2. Bounds outside the two root arcs

We use Wang--Krattenthaler,
[arXiv:2201.12415, Lemmas 9.2 and 9.4](https://arxiv.org/html/2201.12415#S9).
These are bounds for $P_n$ itself and do not restrict the power to
$\{1,2,3\}$. Their bounds for the loss
$D=\log|P_n(\omega r)/P_n(re^{i\theta})|$ are

```math
D\ge-\frac{0.8}{|1-re^{i\theta}|}
 +\frac{r^3(1+r^3)}6\frac{1-r^{n/2}}{1-r^3}
       \left(1-(1+18\rho^2)^{-1/2}\right),              (C3)
```

when $\theta$ is within $\pi/3$ of a multiple of $2\pi/3$,
with $\rho=(\theta-2\pi h/3)(1-r^{3n})/(1-r^3)$; and

```math
D\ge\frac{(r+r^2)(1-r^{3n})}{6(1-r^3)}-5.44          (C4)
```

when $|1-re^{i\theta}|<1/3$. At $r=1$ the ratios have their limits.

Take $r=e^{-\tau/(3n)}$, $0\le\tau\le7$, and choose a fixed,
sufficiently small $\eta>0$. Exclude arcs of half-width
$\eta/(3n)$ around the two nontrivial roots. On the remaining circle,
if $|1-re^{i\theta}|<1/3$, (C4) gives $D\ge c n-C$, uniformly.
Otherwise, the negative term in (C3) is at least $-2.4$.
If the nearest root is nontrivial, the excluded-arc condition gives
$|\rho|\ge c_\eta>0$. If the nearest root is 1, then
$|1-re^{i\theta}|\ge1/3$, together with $r\to1$ uniformly,
keeps its angular distance bounded below, so the same lower bound holds.
Finally $(1-r^{n/2})/(1-r^3)\ge c n$ uniformly on this compact
$\tau$-interval. Hence (C3) again gives $D\ge c n-C$.
Raising the modulus ratio to the fourth power proves an exponentially
small bound outside the two arcs. Constants can depend on $\eta$, but
are independent of $n,\tau,j$.

## 3. The expansion to the precision needed for adjacent coefficients

For the saddle $\tau$ of $j$, put $g=H''(\tau)>0$, and

```math
\mathcal B_{n,j}=
 \frac{\exp\{n(H(\tau)+\alpha(\tau)\tau)\}}
 {3n\sqrt{2\pi n g}}>0.                              (C5)
```

For $j\equiv a\pmod3$, uniformly when $0\le\tau\le7$,

```math
\frac{c_j(n)}{\mathcal B_{n,j}}
=\Psi_a(\tau)+\frac{A_{1,a}(\tau)}n
 +\frac{A_{2,a}(\tau)}{n^2}+O(n^{-3}).                 (C6)
```

The functions $A_{1,a},A_{2,a}$ are real analytic on a neighborhood
of this real interval. In particular they and their first derivatives
are bounded there. The first is

```math
A_{1,a}=b\Psi_a-\frac{\Psi_a''}{2g}
 +\frac{H'''\Psi_a'}{2g^2}
 +\Psi_a\left(\frac{H''''}{8g^2}-\frac{5(H''')^2}{24g^3}\right).  (C7)
```

Here is a justification of the remainder order, including why an
ordinary value estimate suffices. On either arc put
$q=\omega^{\pm1}e^{-(\tau-iy)/(3n)}$. The measure is
$dy/(3n)$, and the coefficient phase is $e^{-ijy/(3n)}$.
Combining the two root weights gives $\Psi_a$, with the common
exponent $n[H(\tau-iy)-H(\tau)-i\alpha y]$.
The linear term vanishes. Its quadratic term is $-ngy^2/2$.
After reducing $\eta$ if necessary, analyticity and the positive
minimum of $g$ give a uniform real-part bound $-cny^2$.

Set $v=\sqrt n\,y$. In the region $|v|\le n^{1/20}$, expand
the exponential and the amplitude from (C1) in $n^{-1/2}$ through
degree 5. Uniform Taylor remainders are bounded by
$C n^{-3}(1+|v|^K)e^{-cv^2}$, for a fixed finite $K$, after
including the Gaussian weight. This follows from the analytic Taylor
remainder and $n|y|^3=|v|^3/\sqrt n=o(1)$ there. All odd-degree
terms are odd in $v$, so their symmetric integrals vanish. The even
terms of degree 0, 2, 4 define $\Psi_a,A_{1,a},A_{2,a}$.
Outside this region, the Gaussian bound is exponentially small in
$n^{1/10}$, and polynomial Taylor terms have the same property.
The exterior-circle estimate in Section 2 contributes
$O(n^{3/2}e^{-cn})$ after (C5). This proves (C6).
The usual Gaussian moments give (C7).

For an explicit construction of the second coefficient, let
$\mathcal S_1 f$ be the right side of (C7) with $b\Psi_a$
removed and $\Psi_a$ replaced by $f$. Define

```math
\mathcal S_2 f=
\sum_{p+e_3+2e_4+3e_5+4e_6=4}
 \frac{(-1)^{D/2}(D-1)!!}{p!\,g^{D/2}} f^{(p)}
 \prod_{r=3}^6\frac{(H^{(r)})^{e_r}}{(r!)^{e_r}e_r!},
\quad D=p+3e_3+4e_4+5e_5+6e_6.
```

The indices in this finite sum are nonnegative integers; $D$ is even.
Then

```math
A_{2,a}=\frac{b^2}{2}\Psi_a+
 d(\omega^{-a}e^h-\omega^ae^{-h})
 +\mathcal S_1(b\Psi_a)+\mathcal S_2\Psi_a.            (C7a)
```

This also makes its analyticity and bounded derivative explicit. The
finite-sum expression is evaluated against exact coefficients
by `saddle_diagnostic.py`; that comparison is diagnostic only.

## 4. One change, including the tiny transition window

Take a fixed closed neighborhood $J\subset(0,6)$ of $\tau_*$
on which $\Psi_2'>0$. For two consecutive residue-2 degrees
$j,j+3$ whose saddles lie in $J$, the saddle equation gives

```math
\tau(j+3)-\tau(j)=-\frac{1}{n^2H''(\tau(j))}+O(n^{-4}).
```

Subtract (C6) for those two degrees. Analyticity of $A_1,A_2$
bounds their differences, while the two separate remainders together
are still $O(n^{-3})$. Thus

```math
\frac{c_{j+3}(n)}{\mathcal B_{n,j+3}}
-\frac{c_j(n)}{\mathcal B_{n,j}}
=-\frac{\Psi_2'(\tau(j))}{n^2H''(\tau(j))}+O(n^{-3})<0          (C8)
```

for all sufficiently large $n$, uniformly in $J$. This avoids
assuming that an error estimate can simply be differentiated.
Outside $J$, (C6) has the strict sign of $\Psi_a$, since its
nonzero leading amplitudes have uniform margins on the remaining
compact intervals. Consequently, in the compact saddle range the
residue-2 sequence changes sign exactly once, with at most one zero.
The normalization in (C5) is positive at every lattice point, so it
preserves coefficient signs. This statement is joined to the endpoint
and initial ranges in LARGE_N_RESULT.md.

## 5. Location, with a linear correction

The same argument yields a last nonnegative dissection index $k_n$
with

```math
k_n=\alpha_*n^2+\beta_*n+O(1),\qquad
\beta_*=-\frac{\Psi_2''(\tau_*)}{2\Psi_2'(\tau_*)}
          +\frac{H'''(\tau_*)}{2H''(\tau_*)}.           (C9)
```

Indeed set $\tau=\tau_*+s/n$ in (C6). Its order-$1/n$ zero
has $s=-A_{1,2}(\tau_*)/\Psi_2'(\tau_*)$. Expanding
$j/(3n^2)=-H'(\tau)$ gives (C9); $j=3k+2$ changes only the
bounded term. The residual normalized error is $O(n^{-2})$, and
(C8) is bounded below in magnitude by a positive multiple of
$n^{-2}$. Therefore the uncertainty in the integer crossing index
is bounded, as asserted.

`certify_limit.py` encloses $\alpha_*$ by exact rational interval
arithmetic:

```math
0.74908009478851072443256<\alpha_*<0.74908009478851072443257.
```

It uses $\tan(\pi/24)=\sqrt6+\sqrt2-\sqrt3-2$, integer square-root
enclosures, Machin's arctangent formula for $\pi$, and the identity

```math
R(\tau)=\frac{\pi^2/9+\operatorname{Li}_2(e^{-3\tau})/3
                      -\operatorname{Li}_2(e^{-\tau})}{\tau}.
```

The logarithm and dilogarithm are enclosed by positive series and
geometric tail bounds. The linear correction is currently a high-precision
diagnostic value, $\beta_*=0.0881422570523141372\ldots$.
The exact expressions define both constants in the theorem.
