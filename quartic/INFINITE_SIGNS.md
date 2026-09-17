# Infinite fourth-power signs and the initial coefficients

17 September 2026. Derivation and code by OpenAI Codex; anonymous hobby
research. No independent review or formal verification.

Put $F_\infty(q)=((q;q)_\infty/(q^3;q^3)_\infty)^4$, and
$b_j=[q^j]F_\infty$. Then

```math
b_{3k}>0,\qquad b_{3k+1}<0,\qquad b_{3k+2}>0\quad(k\ge0).       (I1)
```

We use the explicit remainder formula in Liuquan Wang,
[arXiv:2108.03932v3, Theorem 1.1](https://arxiv.org/pdf/2108.03932).
Its parameters are $t=3,m=4,\mu=1/3,A=1,M=\sqrt2/3$; its
hypothesis $m(t-1)\le24$ holds. For $d=j-1/3>0$, put
$X=4\pi\sqrt d/(3\sqrt3)$. Substitution into that theorem gives

```math
b_j=\frac{2\pi}{3\sqrt3\sqrt d}\,
 2\cos(2\pi j/3+2\pi/9)I_1(X)+E_j,
\qquad |E_j|<4e^{X/2}+10^7.                              (I2)
```

The cosine signs are $+,-,+$ in residues $0,1,2$. Its smallest
absolute value, including the factor 2, is $2\sin(\pi/18)>1/3$.
For example, $\pi>157/50$ and $\sin u\ge u-u^3/6$ prove the
last strict inequality.

Here are deliberately loose bounds for the constants in (I2). The
coefficient of $e^{X/2}$ is
$(\pi^7/24)^{1/4}<4$. The two constant remainder terms in Wang's
formula are bounded by $6e^{11}$ and $54e^8$, respectively.
For the latter, use $e^{-\pi}<1/20$,
$e^{-\pi/3}<3/8$, and $\pi/3<22/21$: the exponent inside
that term, apart from its initial $e^2$, is less than 6.
Finally $6\cdot3^{11}+54\cdot3^8<10^7$.

For $X\ge1$, the real integral representation
([DLMF 10.32.2](https://dlmf.nist.gov/10.32#E2))

```math
I_1(X)=\frac X\pi\int_{-1}^{1} e^{Xt}\sqrt{1-t^2}\,dt
```

restricted to $[1-X^{-1},1-(2X)^{-1}]$ gives
$I_1(X)\ge e^{X-1}/(2\pi\sqrt{2X})>e^X/(25\sqrt X)$.
Consequently the absolute main term in (I2) exceeds
$e^X/(200d^{3/4})$. We used
$2<X/\sqrt d<5/2$, $\sqrt{5/2}<8/5$, and
$2\pi/(3\sqrt3)>1$.

For $d\ge300$, the error divided by this lower bound is at most

```math
800d e^{-\sqrt d}+2\cdot10^9d e^{-2\sqrt d}
\le240000e^{-17}+6\cdot10^{11}e^{-34}<1/50.           (I3)
```

Both functions on the first line decrease for $d\ge300$.
The last bound is certified by a positive rational Taylor lower bound
for $e^{17}$. Thus (I1) holds for $j\ge301$.
`certify_infinite.py` checks every coefficient $0\le j\le300$ using
exact integers. These checks and (I2)--(I3) prove (I1) for all indices.

**Source-table caution.** The extracted text of Table 3 in the cited
arXiv version lists the power-4 positive and negative sets inconsistently
with its own formula (1.21), the elementary coefficient $b_2=2$, and
(I2). The argument here uses Theorem 1.1 and the explicitly evaluated
cosines, not that table row. No inference is made about later versions.

For the finite polynomial $F_n=P_n^4$, its first omitted factor has
degree $3n+1$. Hence
$[q^j]F_n=b_j$ for $0\le j\le3n$. This proves all required
initial signs uniformly in $n$, without a finite search over $n$.
