# 第三 Borwein 予想の eventual form に向けた統合証明稿

*This document is written in Japanese.*
## 零剰余類の先行除去、端点評価、内部評価と再検証

2026-09-08 / 統合稿 v0.7 / 未査読・新規性未確定

### 要旨と位置づけ

有限積

```math
B_n(q)=\prod_{\substack{1\le j\le5n\\5\nmid j}}(1-q^j)
      =\sum_{m=0}^{10n^2}c_n(m)q^m
```

に対し、十分大きい $n$ で第三 Borwein 予想の符号を全係数に保証するための証明案を、原ノートから一つに再構成する。端の二つの剰余類では無限積の係数が恒等的に零になる。この事実を係数積分の前に利用し、小さな尾部因子を局所主項だけでなく副円弧・eta 剰余にも残す。最初の尾部区間は Rogers–Ramanujan 級数の非負恒等式で扱い、それ以後は位置を $5n$ ずらした鞍点で評価する。内部は有限積そのものの Fourier 平滑化と有理近似により処理する。

本稿の作業定理は $n\ge31147$ の全係数版である。v0.7はv0.6の証明骨格を保ち、分母別パラメータ、等比級数核、共鳴ギャップの丸め、中間円弧の分割、端点の許容幅を調整して4万から下げた。固定した端点上限 $w_0=0.0013$ では、この整数が採用した端点条件の下限である。証明手法全体の最適閾値を主張するものではない。本文は独立監査前の解析的証明案であり、数値比較の区間検証、有限整数検査、通常多倍長診断を区別する。

全 $n$ の第三 Borwein 予想の解決、新規性の確定、第三者による独立監査の完了は主張しない。原ノートで参照されていたプログラム・証明書は受領していないため、付属コードと出力は本統合作業で再構築したものである。

## 1. 主張、記号、証明の依存関係

### 作業定理 A

整数 $n\ge31147$、$0\le m\le10n^2$ に対して

```math
c_n(m)\ge0\quad(5\mid m),\qquad c_n(m)\le0\quad(5\nmid m).
\tag{1.1}
```

同じ $n$ の範囲で、零集合は

```math
\begin{aligned}
Z_n^-&=\{5j+3,5j+4:0\le j<n\}\cup\{7,5n+8,5n+9\},\\
Z_n&=Z_n^-\cup(10n^2-Z_n^-),\qquad |Z_n|=4n+6.
\end{aligned}
\tag{1.2}
```

以下ではこの主張を導く証明を提示する。未査読の作業定理という位置づけは全節に適用する。

| 段階 | 内容 | 主な依存先 |
|---|---|---|
| §2–3 | 無限積の5分解、非負核、最初の尾部公式 | 古典恒等式、形式的冪級数 |
| §4–5 | 端の明示的相対誤差 | eta 変換、尾部展開、全角度ギャップ |
| §6 | $n\ge31147$ の端領域 | §3、§5 |
| §7–9 | 固定内部区間の絶対誤差 | 位相下界、有限積の局在化、Gaussian 積分 |
| §10 | 全係数と零集合の接続 | §6、§9、相反対称性 |

原ノートの非明示的な eventual theorem を前提にはしない。小素数全体に関する先行ノートの定理にも依存せず、必要な $p=5$ の内部評価は本稿に記す。

$D=10n^2$、$\zeta=e^{2\pi i/5}$ とする。因子数は $4n$ なので

```math
c_n(m)=c_n(D-m).
\tag{1.3}
```

反転は剰余類 $a$ を $-a\bmod5$ に移すが、所望符号は保存する。

基本関数は

```math
\begin{aligned}
h(z)&=\frac{\operatorname{Li}_2(e^{-z})-\pi^2/6}{z},\\
R(z)&=h(5z)-h(z),\qquad
C=\frac{2\pi^2}{15},\quad A=\frac C5.
\end{aligned}
\tag{1.4}
```

$h$ は右半平面から枝を取る。実 $\tau\ge0$ では

```math
R(\tau)=\int_0^1\log\left(\sum_{j=0}^4e^{-j\tau x}\right)dx,\quad
R(0)=\log5,\quad R'(0)=-1,\quad R''(0)=\frac23.
\tag{1.5}
```

重み $\Pr(J=j)\propto e^{-j\tau x}$ を使うと

```math
R''(\tau)=\int_0^1x^2\operatorname{Var}(J)\,dx>0.
\tag{1.6}
```

従って $-R'$ は1から0へ厳密に減少する。$0<k<5n^2$ なら

```math
-R'(\tau)=\frac{k}{5n^2}
\tag{1.7}
```

は唯一の正の解を持つ。$k=5n^2$ は $\tau=0$ で扱う。

## 2. 無限積の5分解と零剰余類

$(a;Q)_j=\prod_{r=0}^{j-1}(1-aQ^r)$ と定める。

```math
G(q)=\frac{(q;q)_\infty}{(q^5;q^5)_\infty},\qquad
T_n(q)=\prod_{\substack{j>5n\\5\nmid j}}(1-q^j)^{-1},\qquad B_n=GT_n.
\tag{2.1}
```

Rogers–Ramanujan 級数を

```math
R_t(Q)=\sum_{j\ge0}\frac{Q^{j^2+(t-1)j}}{(Q;Q)_j},\qquad
g=R_1,\quad h_R=R_2
\tag{2.2}
```

とする。古典的な和・積恒等式 [RR] により

```math
g(Q)=\frac1{(Q,Q^4;Q^5)_\infty},\qquad
h_R(Q)=\frac1{(Q^2,Q^3;Q^5)_\infty}.
```

古典的な5分解 [W21, Proposition 3.6 の証明] は

```math
G(q)=g(q^5)^2-qg(q^5)h_R(q^5)-q^2h_R(q^5)^2.
\tag{2.3}
```

同文献の連分数表示からは、$Q=q^5$、連分数の積比を $h_R(Q)/g(Q)$ と置き、$(q^{25};q^{25})_\infty/(q^5;q^5)_\infty=g(Q)h_R(Q)$ を掛けることでこの形が出る。

右辺の三つの $Q$ 級数は係数非負であり、

```math
[q^{5j+3}]G=[q^{5j+4}]G=0.
\tag{2.4}
```

従って $m\equiv3,4\bmod5$ なら、厳密に

```math
c_n(m)=[q^m]\{G(q)(T_n(q)-1)\}.
\tag{2.5}
```

(2.5) は円周全体の係数積分に使用する。主円弧の主項だけから $G$ を引く操作ではない。

なお、$g^2,gh_R$ の全係数は正である。$h_R^2$ は $Q^1$ のみ零で、$Q^0$ と $Q^K$（$K\ge2$）は正である。後者は部品2と3で全整数 $K\ge2$ を作れることから従う。従って $[q^7]G=0$ であり、小指数も含めた厳密な負号を無条件に引用しない。

## 3. 非負核と最初の尾部係数の厳密公式

定義から、差の分子の $1-Q^j$ を消し $j=i+1$ と置けば

```math
R_t=R_{t+1}+Q^tR_{t+2}.
\tag{3.1}
```

$W_t=R_{t+1}^2-QR_tR_{t+2}$ と置き、(3.1) を代入すると

```math
W_t=(1-Q)R_{t+1}R_{t+2}+Q^{2t+2}W_{t+2}.
\tag{3.2}
```

反復して残項の最低次数を無限大へ送ることにより、形式的冪級数として

```math
\frac{R_{t+1}^2-QR_tR_{t+2}}{1-Q}
=\sum_{j\ge0}Q^{2j(t+j)}R_{t+2j+1}R_{t+2j+2}.
\tag{3.3}
```

特に $QR_3=g-h_R$ より

```math
\mathcal D(Q):=\frac{gh_R+h_R^2-g^2}{1-Q}
=\sum_{j\ge0}Q^{2j(j+1)}R_{2j+2}R_{2j+3}
=\sum_{K\ge0}d_KQ^K.
\tag{3.4}
```

右辺は非負で、$d_0=1,d_1=0$。$j=0$ 項は $Q^2/(1-Q)$ を含むので $d_K>0$（$K\ge2$）。

尾部の因子から選ぶ個数で展開すると

```math
T_n(q)=1+q^{5n}\frac{q+q^2+q^3+q^4}{1-q^5}+O(q^{10n+2}).
\tag{3.5}
```

(2.3) と掛けたとき、剰余類3、4のどちらにも $(g^2-gh_R-h_R^2)/(1-Q)=-\mathcal D$ が現れる。よって

```math
\boxed{c_n(5n+5K+a)=-d_K\quad(a=3,4,\ 0\le K\le n-1).}
\tag{3.6}
```

対象次数は $10n-1$ 以下なので、(3.5) の剰余は寄与しない。これは全 $n$ の公式であり、$K\ge n$ へは延長していない。

## 4. 端点の主項と明示的誤差

$\tau,w_0,x_0$ を

```math
-R'(\tau)=\frac{k}{5n^2},\qquad
w_0=\frac{\tau}{5n},\qquad x_0=e^{-\tau}
```

で定め、

```math
\tau\ge11/2,\qquad 0<w_0\le0.0013
\tag{4.1}
```

を仮定する。正の主項と誤差予算は

```math
P(n,k)=
\frac{\exp\{nR(\tau)+kw_0-w_0/6\}}
 {5\sqrt{2\pi R''(\tau)}\,n^{3/2}},
\tag{4.2}
```

```math
\mathcal E(\tau,w_0)=15\sqrt{w_0}+15e^{-\tau}+60w_0
+50w_0^{-3}e^{-1/(30w_0)}.
\tag{4.3}
```

### 端点評価 E

弱い二剰余類では $k=m-5n$、強い三剰余類では $k=m$ とすると

```math
\left|\frac{c_n(m)}{-P(n,m-5n)}-1\right|\le\mathcal E
\quad(m\equiv3,4\bmod5),
\tag{4.4}
```

```math
\left|\frac{c_n(m)}{\sigma_aP(n,m)}-1\right|\le\mathcal E
\quad(m\equiv a\bmod5,\ a=0,1,2),
\tag{4.5}
```

```math
(\sigma_0,\sigma_1,\sigma_2)
=\left(\frac{5+\sqrt5}{2},-\sqrt5,-\frac{5-\sqrt5}{2}\right).
```

$\tau$ について誤差は減少し、$w_0$ について増加する。最後の項の単調性は $1/30-3w_0>0$ による。従って

```math
\mathcal E\le \mathcal E(11/2,0.0013)
=0.84662267115518262237\ldots<0.85<1.
\tag{4.6}
```

次節で (4.4)–(4.5) を導く。

## 5. 端点評価の証明

### 5.1. 正係数周期積の全角度ギャップ

$\mathscr P\in\{g^2,gh_R,h_R^2\}$ とする。積の重みは周期5で、非零剰余類の重みはそれぞれ

```math
(2,0,0,2),\quad(1,1,1,1),\quad(0,2,2,0).
```

$\sum b_r=4,\sum rb_r=10$ である。$0<v\le0.0065$ と
$\operatorname{dist}(\theta,2\pi\mathbb Z)\ge3v/4$ のもとで

```math
|\mathscr P(e^{-v+i\theta})|
\le\mathscr P(e^{-v})e^{-1/(5v)}.
\tag{5.1}
```

実際、対数級数の第一項だけで損失は

```math
S=\sum_{k\ge1}b_ke^{-vk}(1-\cos k\theta)
=\frac{N_v(0)}{1-e^{-5v}}
-\Re\frac{N_v(\theta)}{1-e^{-5v+i5\theta}},
\quad N_v(\theta)=\sum_{r=1}^4b_re^{-rv+ir\theta}
```

以上となる。$\theta=2\pi b/5+\epsilon,\ |\epsilon|\le\pi/5$ と書く。

$|\epsilon|\ge2v$ では、分母の平方と
$1-\cos t\ge2t^2/\pi^2$（$|t|\le\pi$）から

```math
S\ge\frac{4-10v}{5v}
\left(1-\frac1{\sqrt{1+16e^{-5v}/\pi^2}}\right)\ge\frac{0.2}{v}.
```

$|\epsilon|\le2v$ では $t=\epsilon/v$、
$C_b=\sum b_re^{2\pi ibr/5}$ と置く。対称な重みなので $C_b$ は実数、
$C_0=4$、$b\ne0\bmod5$ では $C_b\le\sqrt5-1<1.25$。

```math
S=\frac{4-C_b/(1+t^2)}{5v}+E,\qquad |E|<13.
\tag{5.2}
```

ここで $|N_v(\theta)-C_b|\le30v$、
$|1-e^{-5v+i5\epsilon}|\ge5v/(1+5v)$ を使う。
$d=5v-i5\epsilon,\ |d|\le5\sqrt5v$ に対して [Coth] から

```math
\left|\frac1{1-e^{-d}}-\frac1d\right|
\le\frac12+\frac{|d|}{12(1-|d|^2/(4\pi^2))}<0.507.
```

非零角側の分子置換と分母置換の誤差は $6.195+2.028$ 以下、
実軸側の誤差は $2+4(0.507)=4.028$ 以下で、総和は13未満。
$b=0$ なら $|t|\ge3/4$ なので主項は $0.288/v$ 以上、
$b\ne0$ なら $0.55/v$ 以上。$13v\le0.0845$ を引いても
$0.2/v$ 以上となり (5.1) が従う。

この議論は角度をサンプルする検査ではない。スカラー比較は証明書 E01–E06 が扱う。

### 5.2. eta 変換と正の実軸の上界

[Eta] の変換則を $q=\zeta^\ell e^{-w}$ に適用すると

```math
G(\zeta^\ell e^{-w})
=\kappa_\ell e^{A/w-w/6}U_\ell(w),\quad
(\kappa_1,\kappa_2,\kappa_3,\kappa_4)
=(e^{-i\pi/5},1,1,e^{i\pi/5}),
\tag{5.3}
```

```math
U_\ell(w)=
\frac{(\zeta^{-\ell^{-1}}e^{-4\pi^2/(25w)};
        \zeta^{-\ell^{-1}}e^{-4\pi^2/(25w)})_\infty}
     {(e^{-4\pi^2/(5w)};e^{-4\pi^2/(5w)})_\infty}.
\tag{5.4}
```

$\ell^{-1}$ は5を法とする逆元。乗数には Dedekind 和
$s(\ell,5)=(1/5,0,0,-1/5)$ を用いる。
また

```math
G(e^{-w})=\sqrt5 e^{-C/w-w/6}U_0(w),\quad
U_0(w)=\frac{(e^{-4\pi^2/w};e^{-4\pi^2/w})_\infty}
 {(e^{-4\pi^2/(5w)};e^{-4\pi^2/(5w)})_\infty}.
\tag{5.5}
```

$w=w_0+iy,\ |y|\le3w_0/4$ では
$\Re(1/w)\ge16/(25w_0)$。
$|\log(z;z)_\infty|\le |z|/(1-|z|)^2$ と
$(4\pi^2/25)(16/25)>1$ より、(4.1) の範囲で

```math
|U_\ell(w)-1|\le10e^{-1/w_0}\quad(\ell=0,\ldots,4).
\tag{5.6}
```

(2.3) の根のフィルターで三つの $\mathscr P$ を取り出す。
$w=v/5$ として、原始根の各項は $1.001e^{A/w}$ 以下、
$q=1$ の項は $0.001e^{A/w}$ 以下。$q,q^2$ の係数を外す因子は
高々 $e^{2w}$ なので

```math
\mathscr P(e^{-v})\le2e^{C/v}\qquad(0<v\le0.0065).
\tag{5.7}
```

証明書 E07–E10、E36–E37 が数値余裕を確認する。

### 5.3. 尾部展開と四つの円弧の合成

$x=e^{-5nw}$、$|x|=x_0$ として

```math
\Lambda(x)=\frac15\operatorname{Li}_2(x^5)-\operatorname{Li}_2(x),\qquad
L_\ell(x)=\sum_{\substack{j\ge1\\5\nmid j}}\frac{x^j}{j}
\left(\frac{\zeta^{\ell j}}{1-\zeta^{\ell j}}+\frac12\right).
```

正確な尾部対数は

```math
\log T_n(\zeta^\ell e^{-w})
=\sum_{j\ge1}\frac{x^j}{j}
\left\{\frac1{\zeta^{-\ell j}e^{jw}-1}-\frac1{e^{5jw}-1}\right\}.
```

これを展開すると

```math
\log T_n(\zeta^\ell e^{-w})
=\frac{\Lambda(x)}{5w}+L_\ell(x)+E_\ell(w,x),
\qquad |E_\ell(w,x)|\le\frac{5|w|x_0}{1-x_0}.
\tag{5.8}
```

誤差の一様性を確認する。$|\Im t|\le3\Re t/4$ では非自明な5乗根 $\xi$ に対し

```math
\left|\frac{d}{dt}\frac1{\xi^{-1}e^t-1}\right|
=\frac1{2(\cosh\Re t-\cos(\Im t-\arg\xi))}\le4.
```

$\Re t\le1$ は角距離 $2\pi/5-3/4$、$\Re t\ge1$ は $\cosh1-1$ を使う。
また $\Re(t^2)\ge0$ と coth の部分分数展開から

```math
\left|\frac1{e^t-1}-\frac1t+\frac12\right|\le\frac{|t|}{12}.
```

よって各括弧の一次剰余は $5\nmid j$ で $5j|w|$ 未満、
$5\mid j$ で $j|w|/2$ 以下。級数を足して (5.8) を得る。
ここで $x_0/w_0\to0$ は仮定していない。

$C_a=\sum_{\ell=1}^4\zeta^{-a\ell}\kappa_\ell$ とすると

```math
(C_0,C_1,C_2,C_3,C_4)=(\sigma_0,\sigma_1,\sigma_2,0,0).
```

弱い二剰余類について $S_a(x)=\sum\zeta^{-a\ell}\kappa_\ell e^{L_\ell(x)}$ は
$S_a(0)=0,S_a'(0)=-1$ を満たす。導関数の計算は次で確認できる。
$V_a=\sum\zeta^{-a\ell}\kappa_\ell/(1-\zeta^\ell)$ とすると
$V_a-V_{a-1}=C_a$、$\sum_aV_a=0$。
$C_3=C_4=0$ より $V_2=V_3=V_4=-1$ となり、
$S_a'(0)=V_{a-1}+C_a/2=-1$。

$L_\ell$ の $j$ 次係数の絶対値は $1/j$ 以下なので

```math
|S_a(x)+x|\le8x_0^2.
```

(5.8) の指数化による四根合計の誤差は $25|w|x_0$ 以下である。
$(C+\Lambda(x))/(5w)=nR(5nw)$ を使うと

```math
\sum_{\ell=1}^4\zeta^{-a\ell}\kappa_\ell e^{A/w-w/6}
(T_n(\zeta^\ell e^{-w})-1)
=-x e^{nR(5nw)-w/6}(1+\epsilon(w)),
\tag{5.9}
```

```math
|\epsilon(w)|\le8x_0+25|w|.
```

$-1$ の項は $C_a=0$ によって正確に消える。eta 剰余は必ず
$(G-G_{\rm main})(T_n-1)$ として残す。強い三剰余類は
$S_a(0)=\sigma_a\ne0$ を使い、同じ相対誤差上界に収まる。
例えば $|S_a(x)-\sigma_a|/|\sigma_a|\le8x_0$ は
$\min|\sigma_a|=(5-\sqrt5)/2>1$ から従う。
数値比較は E11–E15、E35。

### 5.4. 副円弧・eta 剰余の積分

$r=e^{-w_0}$ とする。尾部の係数非負性から、円周全体で

```math
|T_n(q)-1|\le T_n(r)-1\le\frac{x_0}{w_0}e^{x_0/w_0}.
\tag{5.10}
```

実際、尾部対数を等比級数で評価すると

```math
\log T_n(r)\le
4\frac{x_0}{1-x_0}+\frac4{5w_0}\frac{x_0}{1-x_0}
\le\frac{x_0}{w_0}.
```

五つの5乗根から角距離 $3w_0/4$ 以上の領域では、
$Q=q^5$ に (5.1)、(5.7) を適用して

```math
|G(q)|\le6\exp\left(\frac A{w_0}-\frac{0.04}{w_0}\right).
\tag{5.11}
```

弱い場合は (5.10) の因子を保持するので
$x_0r^{-m}=r^{-k}$（$k=m-5n$）。
$\Lambda(x_0)\ge-x_0/(1-x_0)$ と

```math
0.04-x_0-\frac{x_0}{5(1-x_0)}>\frac1{30}
```

により、(4.2) で割った副円弧寄与は

```math
16w_0^{-3}e^{-1/(30w_0)}
```

以下となる。具体的には §5.5 の $E''(w_0)\le0.67/w_0^3$ より、
前因子は $6\sqrt{2\pi\cdot0.67}\,e^{w_0/6}w_0^{-5/2}$
で、これを $16w_0^{-3}$ 以下に抑える。

$q=1$ の近傍は (5.5)、原始根の eta 剰余は (5.6) を使う。
それぞれ主項との指数差を $0.9/w_0$ 以上にでき、粗い前因子を
$100w_0^{-5/2}$ としても、各々
$w_0^{-3}e^{-1/(30w_0)}$ に吸収できる。
従って除外部分全体には、余裕を取って

```math
50w_0^{-3}e^{-1/(30w_0)}
\tag{5.12}
```

を割り当てる。強い場合は $T_n(r)\le e^{x_0/w_0}$ を使い、
$k=m$ で同じ上界になる。比較は E16–E21。

### 5.5. Gaussian 置換と前因子

```math
F(w)=nR(5nw)=\frac A w+Q_n(w),\quad
Q_n(w)=\frac1{5w}\sum_{j\ge1}\frac{a_je^{-5njw}}{j^2},
```

ここで $a_j=-1$（$5\nmid j$）、$a_j=4$（$5\mid j$）。
$r=2,3$、$|y|\le3w_0/4$ で

```math
w_0^{r+1}|Q_n^{(r)}(w_0+iy)|
\le\frac45r!\sum_{j\ge1}\frac{x_0^j}{j^2}
\sum_{h=0}^r\frac{(\tau j)^h}{h!}.
\tag{5.13}
```

$\tau^he^{-\tau j}$ は $\tau\ge11/2,h\le3$ で減少する。$T_*=11/2$、$x_*=e^{-T_*}$ として

```math
B_2=\frac85(1+T_*+T_*^2/2)\frac{x_*}{1-x_*}<0.142,
```

```math
B_3=\frac{24}5\left[
(1+T_*+T_*^2/2)\frac{x_*}{1-x_*}+\frac{T_*^3}{6}\frac{x_*}{(1-x_*)^2}\right].
```

$E(w)=F(w)+kw$ と置くと、鞍点条件 $E'(w_0)=0$ および

```math
\frac{0.38}{w_0^3}\le E''(w_0)\le\frac{0.67}{w_0^3},\qquad
|E'''(w_0+iy)|\le\frac{2.56}{w_0^4}
\tag{5.14}
```

が従う。$A/w$ の実部を正確に計算し、$Q_n$ の二階項を抑えると

```math
\Re E(w_0+iy)-E(w_0)\le-\frac{0.097y^2}{w_0^3},
\tag{5.15}
```

なぜなら $A/(1+(3/4)^2)-B_2/2>0.097$。

二つの指数を直線補間して

```math
\left|e^{E(w_0+iy)-E(w_0)}-e^{-E''(w_0)y^2/2}\right|
\le\frac{2.56|y|^3}{6w_0^4}e^{-0.097y^2/w_0^3}.
```

これを積分して $J_0=\sqrt{2\pi/E''(w_0)}$ で割ると
$14.9\sqrt{w_0}$ 以下。完全 Gaussian の尾部も含め
$15\sqrt{w_0}$ 以下とできる。

(5.9) の振幅と $e^{-iy/6}$ による誤差は
$8x_0+32w_0$ 以下である。
具体的に $|w|\le5w_0/4$ と $|e^{-iy/6}-1|\le(w_0/8)e^{w_0/8}$ を使うと、
$w_0$ の係数は $31.25+e^{w_0/8}(1+8x_0+31.25w_0)/8<32$。
新しい上限でも $x_0<0.005,w_0<0.002$ を代入すればこの比較が成立する。
絶対積分と $J_0$ の比は $\sqrt{0.67/0.194}<1.86$ 以下なので、
振幅の寄与は $15x_0+60w_0$ 以下。
(5.12) と合わせて (4.3) を得る。比較は E22–E33。

四円弧を同じ $y$ で合成した後、弱い場合の
$-x e^{mw}=-e^{kw}$ により位置が $5n$ ずれる。
この段階で余分な「四根の因子4」をもう一度掛けない。

```math
\frac{e^{E(w_0)-w_0/6}}{2\pi}
\sqrt{\frac{2\pi}{E''(w_0)}}
=\frac{e^{nR(\tau)+kw_0-w_0/6}}
 {5\sqrt{2\pi R''(\tau)}\,n^{3/2}},
```

ここでは $E''(w_0)=25n^3R''(\tau)$ を使用した。
これが (4.2) であり、(4.4)–(4.5) が従う。

## 6. $n\ge31147$ の端領域

```math
\rho=-5R'(11/2)=0.2131128927124522702188745809\ldots.
\tag{6.1}
```

$\Lambda(e^{-\tau})<0$ と
$\frac{d}{d\tau}\Lambda(e^{-\tau})
=\log((1-e^{-5\tau})/(1-e^{-\tau}))>0$ より

```math
0<-R'(\tau)<C/\tau^2,\qquad w_0^2<C/(5k).
\tag{6.2}
```

$k\ge5n,\ n\ge31147$ なら

```math
w_0^2<C/(25n)\le(0.0013)^2,
```

比較 E34 は $C/(25(0.0013)^2)<31147$ を確認する。

$0\le m\le\rho n^2$ を考える。

- $a=0,1,2$：$m\le5n$ は $G$ の係数と一致する。$m>5n$ は $k=m$ として端点評価を適用する。
- $a=3,4$：$m\le5n$ は零。$0<k=m-5n<5n$ は (3.6)。$k\ge5n$ は端点評価を適用する。

どの場合も鞍点が必要なときは $k\le m\le\rho n^2$ なので $\tau\ge11/2$。
従って左端の全係数が所望符号を持つ。右端は (1.3) で従う。
この範囲の零係数は (1.2) に挙げた位置に限られる。

$\rho$ の有限桁小数を上向きに丸めて適用範囲を増やしてはならない。
領域の定義は (6.1) の厳密な関数値である。

## 7. 内部主位相の符号と定量下界

原始根の振幅を

```math
\lambda_\ell(z)=\sum_{j=1}^4\left(\frac j5-\frac12\right)
\{\log(1-\zeta^{\ell j}e^{-z})-\log(1-\zeta^{\ell j})\},\qquad
\mathcal A_\ell(z)=e^{\lambda_\ell(z)}
```

と定義する。枝は $z=0$ から連続に取る。実軸上では共役な $j,5-j$ を組にして
$|\mathcal A_\ell(\tau)|=1$。

```math
\Phi_a(e^{-\tau})=\sum_{\ell=1}^4\zeta^{-a\ell}\mathcal A_\ell(\tau).
\tag{7.1}
```

$x=e^{-\tau}$、$v=(1-x)/(1+x)$、

```math
\alpha=\arctan(v\cot(\pi/5)),\quad
\beta=\arctan(v\cot(2\pi/5)),\quad
U=(\alpha+2\beta)/5,\quad V_0=(2\alpha-\beta)/5
```

とする。共役根を組にし、二つの余弦の和を積に直すと

```math
\Phi_a(x)=4\cos(3\pi a/5+U)\cos(-\pi a/5+V_0).
\tag{7.2}
```

$U,V_0$ は $v\in[0,1]$ で0から $\pi/10$ へ厳密に増加する。
例えば

```math
V_0'\ge v_*=\frac{\sin(2\pi/5)-\cot(2\pi/5)}5>0,\quad
U'\ge u_*=\frac{\sin(\pi/5)\cos(\pi/5)+2\sin(2\pi/5)\cos(2\pi/5)}5>0.
```

角度範囲から $\Phi_0>0$、$\Phi_a<0$（$a\ne0,x>0$）。
弱い二剰余類では

```math
-\Phi_3=4\cos(\pi/5-U)\sin(\pi/10-V_0),\quad
-\Phi_4=4\sin(\pi/10-U)\cos(\pi/5+V_0).
```

$\sin y\ge2y/\pi$、$1-v\ge x$ を使い、
$8\cos(\pi/5)v_*/\pi>1/4$、
$8\cos(3\pi/10)u_*/\pi>1/4$ を得る。
強い三剰余類には一定の正の下界がある。よって

```math
\epsilon_a\Phi_a(x)\ge x/4,\qquad
\epsilon_0=1,\quad\epsilon_1=\cdots=\epsilon_4=-1.
\tag{7.3}
```

数値比較は B01–B05。

## 8. 分母別パラメータによる全角度局在化

$n\ge N=31147$、$0\le\tau\le11/2$、$r=e^{-\tau/(5n)}$ とする。
Dirichlet の近似では全てに共通の $Q=4139$ を使い、既約 $a/b$、$b\le Q$、
$|\theta/(2\pi)-a/b|\le1/(bQ)$ を取る。$t=5n(\theta-2\pi a/b)$。
近似を選んだ後で次のパラメータを使うので、円周の被覆は変わらない。

| 分母 | $\eta$ | $K$ | $L$ | 等比級数核の上界 $\kappa_*$ |
|---|---:|---:|---|---:|
| $5\nmid b$ | $11313/1250000$ | 718 | $\log(2/\eta)$ | $\kappa(Z_1)$ |
| $5\mid b$ | $113/40000$ | 2717 | $\log(2/\eta)$ | $\kappa(Z_5)$ |

この節の $\eta,K,L$ は、その行の値を意味する。
$Z_1=(11/2)K/N+10\pi K/Q$、$Z_5=(11/2)K/N+2\pi K/Q$ とし、
§8.3の関数 $\kappa$ を区間演算で直接使用する。丸めた固定定数を誤差予算へ代入しない。
切断尾部は $\delta=10e^{-\eta K}/(\eta K)$、共鳴和の誤差は $4\kappa_*L/n$ とする。

### 8.1. 平滑化した非5倍数分母との競合

$h(z)=\int_0^1\log(1-e^{-zx})dx$、$M(\tau)=\sup_t\Re h(\tau-it)$ とする。
右半平面の最大値原理、原点での $\Re h\to-\infty$、無限遠での $h\to0$ から
$M$ は非負かつ非増加。$f(y)=\log|2\sin(y/2)|$ と置くと、単位円上の平均は
$t^{-1}\int_0^t f(y)dy$。後続周期は最初の周期の平均を0側へ縮める。
$H(t)=tf(t)+\sum_{k\ge1}\sin(kt)/k^2$ の $(\pi,2\pi)$ 内の唯一の零点が最大点である。
1000項と尾部 $1/1000$ で $H(4.97)>0$、$f(4.97)<1/5$ を囲うので $M(0)<1/5$。

平滑化増分では分子の半径も保持する。$0\le\rho\le1$ なら

```math
|1-\rho e^{iv}|\ge\rho|1-e^{iv}|,
\quad
\frac{|1-\rho e^{-\delta+iv}|}{|1-\rho e^{iv}|}
\le1+\frac{\rho(1-e^{-\delta})}{|1-\rho e^{iv}|}
\le1+\frac{\delta}{2|\sin(v/2)|}.
\tag{8.19}
```

最初の不等式の平方差は $(1-\rho)(1+\rho-2\rho\cos v)\ge0$。
特異点は零測度で、対数上界は局所可積分である。
$\mathcal L(u)=(1+u)\log(1+u)-u\log u$ とすると、
$|\sin(v/2)|\ge|v|/\pi$（$|v|\le\pi$）を積分して周期平均は $\mathcal L(\delta/2)$ 以下。
原点から長さ $s\ge1$ の平均は周期平均の $2\pi$ 倍以下、$s\ge2$ なら $\pi$ 倍以下。
一周期未満では周期全体で上から抑え、それ以上では完全周期と残りに分ければ従う。
$\mathcal L$ の凹性から $\mathcal L(bu)/b$ は $b\ge1$ で非増加。

実半径を保持する別の上界は、$u>0$ に対して

```math
\log|1-e^{-u-\delta+iv}|-\log|1-e^{-u+iv}|
\le\log(1+\delta/u).
\tag{8.16}
```

半径変化に関する微分の上界 $1/(e^u-1)\le1/u$ を積分する。
$u=b\tau x,\delta=b\eta$ として積分すると増分は $\mathcal L(\eta/\tau)$ 以下。

以下の共鳴主項を $\mathcal H_b$ と書く。
まず $\tau\le1/2$ では $R(\tau)\ge R(1/2)>1.189$。
$b=1, |t|\le2$ では、全半径に対し $|1-\rho e^{itx}|\le\max(1,2|\sin(tx/2)|)$ より

```math
\mathcal H_1\le4(1-\pi/6)\log(2\sin1).
```

実際 $|t|x\le\pi/3$ では対数は非正、残りの長さは高々 $1-\pi/6$。
$|t|\ge2$ では $\mathcal H_1\le0.8+4\pi\mathcal L(\eta/2)$。
いずれも $R-\mathcal H_1>0.025$。
$b\ge2$ では角度長 $b|t|<1$ なら各因子の絶対値は1以下、
それ以外では $\mathcal H_b\le0.4+4\pi\mathcal L(\eta)$ なので $R-\mathcal H_b>0.1$。

$\tau\ge1/2$ では $A_0=\pi^2/6$、$B(\tau)=\operatorname{Li}_2(e^{-\tau})$ として
$M(\tau)\le B(\tau)^2/(4A_0\tau)$、$R(\tau)\ge(C-B(\tau))/\tau$。
前者は $|\operatorname{Li}_2(e^{-\tau+it})|\le B(\tau)$ と二次式の最大化から従う。
$[1/2,11/2]$ を幅1/40の200区間 $[u,v]$ で覆い、

```math
\frac{C-B(u)-B(u)^2/A_0}{v}-4\mathcal L(\eta/u)>0.025,
\quad
\frac{C-B(u)-B(u)^2/(2A_0)}{v}-2\mathcal L(\eta/u)>0.1
\tag{8.20}
```

を囲った。前者が $b=1$、後者が $b\ge2$ に適用できる。
$B$ は100項と正の幾何級数尾部で上から囲う。
$R(1/2)>1.189$ には $\operatorname{Li}_2(e^{-5/2})$ の正の5項も使用する。
全域で $R>0.238$。これらは L00a–L11 が検査する。

### 8.2. 周波数ごとの非共鳴誤差

$|1-ue^{iv}|\le e^{\eta/2}|1-ue^{-\eta+iv}|$ より

```math
\frac1n\log|B_n(re^{i\theta})|
\le2\eta-\frac1n\Re\sum_{l\ge1}\frac{e^{-\eta l}}l
\left(\sum_{j=1}^{5n}e^{-l\tau j/(5n)+ilj\theta}
-\sum_{j=1}^{n}e^{-l\tau j/n+i5lj\theta}\right).
\tag{8.1}
```

$d=1,5$、$g=\gcd(b,d)$、$B=b/g$ と置く。以下の非共鳴和は必ず **$b\nmid dl$ の部分集合**で取る。
$r_{d,l}=\operatorname{dist}(dla,b\mathbb Z)$ とすれば

```math
\left\|\frac{dl\theta}{2\pi}\right\|\ge\frac{r_{d,l}-dl/Q}{b},
\quad
\left|\sum e^{-\alpha j+idlj\theta}\right|\le\frac{b}{2(r_{d,l}-dl/Q)}.
\tag{8.17}
```

後者は等比級数と単調な正の重みに対するAbel変換による。
距離を $r_{d,l}=g r$ と書くと、$\epsilon=dK/(gQ)<1$ のもとで

```math
\frac1{r-dl/(gQ)}\le\frac1r+\frac{d l}{gQ(1-\epsilon)r^2}.
\tag{8.21}
```

各長さ $B$ のブロックで非零距離 $r$ は高々二回現れる。
減少する正の重みを並べ替える際、欠けた共鳴項には0を補えば上界が保たれる。

```math
S_0=\sum_{l=1}^K\frac{e^{-\eta l}}{l\lceil l/2\rceil},\quad
T(B)=S_0+\frac{2(1+\log K)e^{-\eta B}}{B(1-e^{-\eta B})},\quad
U(B)=\frac{\pi^2}{3(1-e^{-\eta B})}.
\tag{8.22}
```

最初のブロックは $S_0$ 以下。第 $j\ge1$ ブロックの残りは
$e^{-\eta jB}/(jB)$ と $\sum_{i\le K}1/\lceil i/2\rceil\le2(1+\log K)$ で抑える。
二乗逆数は一ブロックで $2\sum_{r\ge1}r^{-2}=\pi^2/3$ 以下。
従って $d$ 側の正規化非共鳴誤差は

```math
\frac{b}{2gn}\left[T(B)+\frac{d}{gQ(1-\epsilon)}U(B)\right]
\tag{8.23}
```

以下。$T,U$ は $B$ の増加とともに減少する。
$5\nmid b,b>2K$ なら、同じ距離を持つ二つの周波数は $l_1\equiv\pm l_2\pmod b$。
異なる二点では $l_1+l_2=b>2K$ が必要で不可能。
この場合は各距離高々一回なので、(8.23)の $T,U$ をそれぞれ
$S_1=\sum_{l\le K}e^{-\eta l}/l^2$、$\pi^2/6$ に置き換えられる。
これが $dK/Q$ を全周波数へ掛ける損失と、大分母の二重出現による損失を除く。

### 8.3. 共鳴和、全分母の分類

$A=l(\tau-it)$、$M=n$ または $5n$ に対して

```math
\sum_{j=1}^M e^{-Aj/M}-M\int_0^1e^{-Ax}dx
=(1-e^{-A})\left(\frac1{e^{A/M}-1}-\frac M A\right).
\tag{8.12}
```

$A=0$ は連続延長。coth の部分分数展開から $|z|<2\pi$ で

```math
\left|\frac1{e^z-1}-\frac1z\right|
\le\kappa(|z|)=\frac12+\frac{|z|}{12(1-|z|^2/(4\pi^2))}.
\tag{8.13}
```

$|A/M|\le(11/2)K/N+10\pi K/(bQ)$。
非5倍数では $b\ge1$、5倍数では $b\ge5$ を使うと、表の $\kappa_*$ が成立する。
$|1-e^{-A}|\le2$ なので二つの共鳴和の誤差は $4\kappa_*L/n$ 以下。
有限和の尾部係数は4、連続主項の尾部係数は6以下なので、尾部合計は $\delta$ 以下。
共鳴主項は

```math
\mathcal H_b=\begin{cases}
\displaystyle\frac4b\int_0^1\log|1-e^{-b\eta-b(\tau-it)x}|dx,&5\nmid b,\\
\displaystyle\frac1j\int_0^1\log\left|\sum_{r=0}^4e^{-rj(\eta+(\tau-it)x)}\right|dx,&b=5j.
\end{cases}
\tag{8.3}
```

$b=5j,j\ge2$ では $\mathcal H_b\le R(j\tau)/j\le R(\tau)/2$。
一般には非5倍数で $4\log2/b$、5倍数で $5\log5/b$ が上界となる。
粗い非共鳴上界 $\gamma_b bL/n$ の係数は、非5倍数で
$[2(1-K/Q)]^{-1}+[2(1-5K/Q)]^{-1}$、5倍数で $0.6/(1-K/Q)$。
$b=5$ は $d=5$ 側が全て共鳴するので $\gamma_5=[2(1-K/Q)]^{-1}$。
$b=1$ は非共鳴項がない。

| 分母 | 主項とのギャップまたは主項上界 | 非共鳴の評価 |
|---|---|---|
| $b=1$ | ギャップ0.025 | 0 |
| $2\le b\le100,\ 5\nmid b$ | ギャップ0.1 | 粗い上界、$b\le100$ |
| $100<b\le1436,\ 5\nmid b$ | 主項 $4\log2/101$ 以下 | (8.23)、$b\le1436,B\ge101$ |
| $1436<b\le Q,\ 5\nmid b$ | 主項 $4\log2/1437$ 以下 | 一回出現の $S_1,\pi^2/6$ |
| $5<b\le300,\ 5\mid b$ | ギャップ $0.238/2$ | 粗い上界、$b\le300$ |
| $300<b\le Q,\ 5\mid b$ | 主項 $5\log5/301$ 以下 | (8.23)、$B_1\ge301,B_5\ge61$ |
| $b=5$、主円弧外 | 次項のギャップ0.00992 | $\gamma_5 5L/n$ |

各行で $2\eta$、共鳴誤差、非共鳴誤差、尾部を引いても $1/1000$ より大きいギャップが残る。
閾値 $N$ で比較し、それ以後は $1/n$ の項と等比級数の引数上界が減少することを使う。
全区間比較は [localization_certificate.json](results/localization_certificate.json) のL12–L19。

### 8.3a. 差ごとの合成による分母5のギャップ

$p_j(z)=e^{-jz}/\sum_{r=0}^4e^{-rz}$、$z=\eta+\tau x$ とする。
個々の積を最悪値へ落とす前に、差 $d$ ごとに

```math
W_d(z)=\sum_{j=0}^{4-d}p_j(z)p_{j+d}(z)
=\frac{\sum_{j=0}^{4-d}y^{2j+d}}{(1+y+y^2+y^3+y^4)^2},\quad y=e^{-z}
\tag{8.24}
```

とまとめる。$S(y)=1+y+\cdots+y^4$ とすると $dW_d/dy=(1-y)Q_d(y)/S(y)^3$。
昇べき順の係数は

| $d$ | $Q_d$ の係数 |
|---|---|
| 1 | $(1,0,0,-4,-7,-7,-4,0,0,1)$ |
| 2 | $(0,2,2,4,2,2,4,2,2)$ |
| 3 | $(0,0,3,4,8,8,4,3)$ |
| 4 | $(0,0,0,4,6,6,4)$ |

$d=2,3,4$ は $y$ で増加、従って $z$ で減少。
$d=1$ では $(1+t)^9Q_1(t/(1+t))$ の係数が
$(1,9,36,80,95,24,-98,-148,-90,-20)$。
符号変化が一回なのでDescartesの符号法則と $Q_1(0)>0,Q_1(1)<0$ により、
$(0,1)$ 内に零点は一つだけ。従って $W_1(z)$ は一山型で、区間最小は端点にある。
これら多項式恒等式は整数演算でも確認する。

$m=400$、$g_{d,i}=\min(W_d(0),W_d(\eta+(11/2)i/m))$、$g_{d,m+1}=0$ とする。
これは単調な下側包絡線である。各特性関数について

```math
-\log\left|\sum_jp_je^{ijtx}\right|
\ge\sum_{d=1}^4W_d(z)(1-\cos(dtx)).
```

階段関数の積分は、$b_i=i/m$ として
$\sum_i(g_{d,i}-g_{d,i+1})b_i(1-\operatorname{sinc}(dtb_i))$ 以上。
$0\le a\le2, |v|\ge a$ で

```math
1-\operatorname{sinc}v\ge a^2/6-a^4/120
\tag{8.8}
```

を使う。最初の半周期ではsincの単調性とTaylor下界、$|v|\ge\pi$ では
$1-1/\pi>8/15$ から従う。
$|t|\ge6/5$、$a_{d,i}=\min(2,(6/5)db_i)$ とすると

```math
\Gamma=\sum_{d=1}^4\sum_{i=1}^{400}(g_{d,i}-g_{d,i+1})b_i
\left(a_{d,i}^2/6-a_{d,i}^4/120\right)>0.00992.
\tag{8.9}
```

区間下端は0.00992より大きい（全桁は証明書に記録）。
半径の平滑化で正係数和の実軸上の値は増えないので $\mathcal H_5\le R-\Gamma$。
主円弧を $|\theta-2\pi\ell/5|\le(6/5)/(5n)$ と取れば、その外で $b=5$ の場合も処理できる。
従って円周の主円弧外全体で

```math
\frac1n\log|B_n(re^{i\theta})|\le R(\tau)-1/1000.
\tag{8.4}
```

全角度は解析的に被覆した。有限角度サンプルによる外挿は使用しない。

## 9. 四根の先行合成と内部積分

$0\le\tau\le11/2$、$z=\tau-it$ とする。小箱を $|t|\le h=2/5$、大箱を $|t|\le6/5$ とする。
非自明な5乗根 $\xi$ と $0\le x\le1$ に対し、$|1-\xi e^{-zx}|$ は小箱で $d=0.75$ より大きく、大箱で0.05より大きい。
角距離 $2\pi/5-h$ と $2\pi/5-6/5$ を使う。

以下では $A_\ell=\mathcal A_\ell=e^{\lambda_\ell}$ と略記する。

### 9.1. Euler–Maclaurin の一次補正と二次剰余

$f\in C^2[0,1]$ に対して

```math
\left|\sum_{r=0}^{n-1}f((r+\alpha)/n)-n\int_0^1f
-B_1(\alpha)(f(1)-f(0))\right|\le\frac1{8n}\int_0^1|f''|.
\tag{9.16}
```

一小区間のPeano核は $(\alpha-u)_+-(1-u)^2/2-B_1(\alpha)(1-u)$ で、絶対値は1/8以下。
各小区間をスケールし、差を望遠鏡和で足せば従う。
$C^3$ の場合、$B_2(\alpha)\Delta f'/(2n)$ も取り出すと

```math
\sum_{r=0}^{n-1}f((r+\alpha)/n)
=n\int_0^1 f+B_1(\alpha)\Delta f+\frac{B_2(\alpha)}{2n}\Delta f'+\varepsilon,
\qquad |\varepsilon|\le\frac1{n^2}\int_0^1|f'''|.
\tag{9.24}
```

三次の核は
$(\alpha-u)_+^2/2-(1-u)^3/6-B_1(\alpha)(1-u)^2/2-B_2(\alpha)(1-u)/2$。
$|B_1|\le1/2,|B_2|\le1/6$ より、その絶対値は $1/2+1/6+1/4+1/12=1$ 以下。
核は二次以下を消す。小区間ごとに足すと(9.24)になる。

$f_j(x)=\log(1-\zeta^{\ell j}e^{-zx})$、$\alpha=j/5$ に適用する。
主積分は $nR(z)$、一次の端点差は§7の $\lambda_\ell(z)$ である。
$b_j=B_2(j/5)$ として

```math
B_n(\zeta^\ell e^{-z/(5n)})
=e^{nR(z)}A_\ell(z)e^{D_\ell(z)/n+\varepsilon_{\ell,n}(z)},
\qquad |\varepsilon_{\ell,n}|\le3400/n^2,
\tag{9.25}
```

```math
D_\ell(z)=\frac z2\sum_{j=1}^4b_j
\left(\frac{\zeta^{\ell j}e^{-z}}{1-\zeta^{\ell j}e^{-z}}
-\frac{\zeta^{\ell j}}{1-\zeta^{\ell j}}\right).
```

小箱で $|z|<5.6$、$|f_j'''(x)|\le2|z|^3e^{-\tau x}/d^3$ なので、
四つを足した積分は $8(5.6)^3/d^3<3400$。
大箱では(9.16)だけを使い、対数剰余は $7200/n$ 以下とする。

### 9.2. 小さい因子を保つ合成

$x=e^{-\tau}$、$M=1.54$ と置く。実軸では $|A_\ell(\tau)|=1$。
小箱では $|\lambda_\ell'|\le0.8x/d$ より $|A_\ell(z)|\le e^{0.4(0.8)/d}<M$。
大箱では、振幅を積表示で直接評価する。正・負の指数の絶対値和はそれぞれ0.4であり、
$1<|1-\xi|\le2$、$0.05<|1-\xi e^{-z}|\le2$ より
$|A_\ell(z)|\le80^{2/5}<6$。

```math
a_1=\frac{0.8Mx}{d},\qquad
a_2=\frac M2\left(\frac{0.8x}{d^2}+\left(\frac{0.8x}{d}\right)^2\right).
\tag{9.18}
```

従って $A_\ell(\tau-it)$ の実軸からの変化は $a_1|t|$ 以下、一次を除いた剰余は $a_2t^2$ 以下。

```math
\Psi_a(z)=\sum_{\ell=1}^4\zeta^{-a\ell}A_\ell(z),\qquad
\Theta_a(z)=\sum_{\ell=1}^4\zeta^{-a\ell}A_\ell(z)D_\ell(z).
\tag{9.26}
```

$\Psi_a(\tau)=\Phi_a(e^{-\tau})$。弱い剰余類 $a=3,4$ では、$x=0$ での定数項が零なので
$\Psi_a(z)=e^{-z}G_a(e^{-z})$ と因子分解できる。
各振幅は $|e^{-z}|<1$ で解析的で、使用する小箱へ解析接続できる。
以下の誤差評価は因子分解だけに依存せず、実軸の正確な $\Phi_a$ と導関数の $x$ 因子を使う。

$b_1=b_4=1/150$、$b_2=b_3=-11/150$。
共役の二項について $u/(1-u)+u^{-1}/(1-u^{-1})=-1$ より

```math
D_\ell(z)=-z/30+ze^{-z}C_\ell(e^{-z}),\qquad
C_\ell(y)=\frac12\sum_{j=1}^4b_j\frac{\zeta^{\ell j}}{1-\zeta^{\ell j}y}.
\tag{9.27}
```

従って $\Theta_a=-z\Psi_a/30+ze^{-z}\sum\zeta^{-a\ell}A_\ell C_\ell$。
一次EM補正でも弱い剰余類の小さい因子を保持できる。
$Z=(\tau^2+h^2)^{1/2}$、$p_a=|\Phi_a(e^{-\tau})|$ として

```math
D_* =Z\left(\frac1{30}+\frac{2x}{25d}\right),\qquad
|\Theta_a|\le\frac Z{30}(p_a+4ha_1)+\frac{8MZx}{25d}=:T_a.
\tag{9.28}
```

$\sum|b_j|=4/25$ を使った。
指数化の残りは各根で

```math
\left|e^{D_\ell/n+\varepsilon_{\ell,n}}-1-D_\ell/n\right|
\le\frac{(3400+D_*^2/2)e^{D_*/N+3400/N^2}}{n^2}.
\tag{9.29}
```

この二次剰余には $x$ の因子を主張しない。明示的な $n^{-2}$ 上界で扱う。

### 9.3. Gaussian 置換と対称積分

実分布 $\Pr(J=j)\propto e^{-j\tau x}$ の分散を $v_x$ とし、
$V=\int_0^1x^2v_xdx$、$W_r=\int_0^1x^rv_xdx$（$r=3,4$）と置く。
$V\le4/3$。§9.4の直接の区間積分は全域で $V>1/80$ を与える。

$Y=J-\mathbb EJ$ とすると $|Y|\le4$、分散は $v_x\le4$。
特性関数の絶対値は平行移動で変わらないので、$J-2\in[-2,2]$ を使えば
$|\mathbb Ee^{i\omega Y}|=|\mathbb Ee^{i\omega(J-2)}|\ge\cos(2|\omega|)\ge\cos0.8>0.69$（$|\omega|\le0.4$）。
$|\mathbb EYe^{i\omega Y}|\le0.4v_x$、二・三・四次のモーメント上界は $v_x,4v_x,16v_x$。
対数の三・四階微分に代入すると係数は

```math
4/.69+4.8/.69^2+2.048/.69^3<23,
\quad16/.69+37.6/.69^2+30.72/.69^3+9.8304/.69^4<240.
```

従って $|R'''|\le23W_3$、$|R^{(4)}|\le240W_4$。
$1-\cos y\ge(1/2-1.6^2/24)y^2$（$|y|\le1.6$）より
$\Re F_n\le-\beta nVt^2$、$\beta=0.39$。

鞍点条件 $-R'(\tau)=m/(5n^2)$ のもとで

```math
P_{\rm bulk}=\frac{e^{nR(\tau)+m\tau/(5n)}}{5\sqrt{2\pi V}\,n^{3/2}},\quad
F_n(t)=n(R(\tau-it)-R(\tau)+iR'(\tau)t),\quad G_0=-nVt^2/2.
```

$D=F_n-G_0$、$C_3=inR'''(\tau)t^3/6$ とすると
$|D|\le(23/6)nW_3|t|^3$、$|D-C_3|\le10nW_4t^4$。
指数間の線分を積分する二次剰余公式から

```math
|e^{F_n}-e^{G_0}-C_3e^{G_0}|
\le\left(10nW_4t^4+\frac{529}{72}n^2W_3^2t^6\right)e^{-\beta nVt^2}.
\tag{9.30}
```

$A(t)=\Psi_a(\tau-it)$、$A_0=\Phi_a$、$A_1=A'(0)$ と書くと

```math
\begin{aligned}
A(t)e^{F_n}-A_0e^{G_0}={}&A_0(e^{F_n}-e^{G_0}-C_3e^{G_0})\\
&+(A(t)-A_0)(e^{F_n}-e^{G_0})
+(A(t)-A_0-A_1t)e^{G_0}+(A_0C_3+A_1t)e^{G_0}.
\end{aligned}
\tag{9.13}
```

最後の項は対称区間で正確に消える。他の項には変化 $4a_1|t|$、二次剰余 $4a_2t^2$ を使う。
実軸全体の偶数Gaussianモーメントを積分し、$\sqrt{2\pi/(nV)}$ で割ると、誤差係数は

```math
H_a=\frac1{\sqrt{2\beta}}\left[
p_a\left(\frac{240W_4}{32\beta^2V^2}+\frac{2645W_3^2}{192\beta^3V^3}\right)
+4\left(\frac{23a_1W_3}{8\beta^2V^2}+\frac{a_2}{2\beta V}\right)\right].
\tag{9.31}
```

最初の括弧に掛かるのは4ではなく、合成した主項の $p_a$ である。
EM一次補正を取り出した積分の寄与は $T_a/(n\sqrt{2\beta})$ 以下。
残るEM剰余の係数を

```math
U=\frac{4M}{\sqrt{2\beta}}(3400+D_*^2/2)e^{D_*/N+3400/N^2}
\tag{9.32}
```

とすれば、その誤差は $U/n^2$ 以下となる。

中間角度 $[0.4,1.2]$ を幅1/20の16区間 $[u,v]$ で覆う。
各区間で $c=\tfrac12(\sin(2v)/(2v))^2$ とすれば、$1-\cos y\ge cy^2$（$|y|\le4v$）。
全16区間で $cu^2>0.052$ を区間検証する。
profileの分散下端から $0.052V>0.00074$、$Vh^2/2>0.001$ を直接検証し、
更新した前因子34を使うと

| 部分 | 正規化した上界 |
|---|---|
| 副円弧 | $15n^{3/2}e^{-n/1000}$ |
| 主円弧の中間角度 | $34\sqrt n\,e^{-37n/50000}$ |
| 完全Gaussianへの延長 | $8e^{-n/1000}$ |

それぞれ $n\ge N$ で $1/n$ 未満。閾値で対数を比較し、導関数 $c-\alpha/n>0$ で以後を保証する。
したがって

```math
\left|\frac{c_n(m)}{P_{\rm bulk}}-\Phi_{m\bmod5}(e^{-\tau})\right|
\le\frac{H_a+T_a/\sqrt{2\beta}+3}{n}+\frac U{n^2}.
\tag{9.33}
```

Cauchy積分の変数は $\theta=2\pi\ell/5+t/(5n)$ であり、
$\sqrt{2\pi/(nV)}/(2\pi\cdot5n)$ が $P_{\rm bulk}$ の前因子を与える。

### 9.4. 全パラメータの区間証明

$\tau_i=(11/2)i/128$、$x_j=j/200$ とし、全長方形上で

```math
v_x=\sum_{j<k}(k-j)^2\frac{e^{-(j+k)\tau x}}{(\sum_{r=0}^4e^{-r\tau x})^2}
```

の値域を区間演算で囲う。これに単項式の正確な積分重みを掛けて足し、$V,W_3,W_4$ を囲う。
位相は(7.2)を同じ区間全体で評価する。各剰余類の位相下端を別々に保持し、
合成後の $p_a$ と誤差の式を区間演算で組み合わせる。

全128区間で、五つの位相の正性、$V>1/80$、五つの誤差・位相比が0.98未満であることを確認した。
計1408比較。分散の下端は約0.01423554、比の最大上端は約0.37280682である。
計算はサンプルの外挿ではなく、全25600長方形の包含評価。
[certify_profiles.py](verification/certify_profiles.py) と [profile_certificate.json](results/profile_certificate.json) に全区間を記録する。
全 $n\ge N$ では右辺の $1/n,1/n^2$ が減少するので、同じ比較が成立する。

## 10. 全係数と零集合への接続

$n\ge31147$ とする。$0\le m\le\rho n^2$、$\rho=-5R'(11/2)$ は§6が覆う。
$\rho n^2\le m\le5n^2$ では未シフト鞍点の解が $0\le\tau\le11/2$ にあり、
§9.4の一様比較により(9.33)の誤差が符号付き主項より小さい。
したがって内部の全係数は厳密に所望符号を持つ。右半分は相反対称性で従う。

零係数は、無限積との一致による $\{5j+3,5j+4:0\le j<n\}\cup\{7\}$、
最初の尾部公式の $d_1=0$ による $5n+8,5n+9$、およびそれらの反転のみである。
他の係数は端点・内部の厳密符号を持つので(1.2)も従う。

この解析的証明案が正しいことを前提としても、$257\le n<31147$ は未処理である。
全$n$版の証明と、解析的推論全体の独立監査は完了していない。

## 11. 再構築したコードと実行結果

### 11.1. 区間証明書

付属コードは mpmath 1.3.0 の区間演算、60桁で次を検査した。

- certify_scalars.py：端点37個、位相5個、内部40個の計82比較。
- certify_localization.py：全分母の分類と半径の全200区間を含む421比較。
- certify_profiles.py：全128区間・1408比較。
- certify_resonant_gap.py：差ごとの多項式恒等式と単調性の根拠を整数で確認し、400区間の一様ギャップを囲う。

十進入力は文字列または整数比で与え、差の区間下端が厳密に正であることを要求する。
通常の浮動小数点への変換を判定に使わない。v0.6で旧版の未使用Bチェックを除き、v0.7では中間円弧の16区間をG16_0–G16_15で検査する。
分散の下端から中間円弧とGaussian尾部の減衰率を直接導く比較はG17–G18、
前因子のN依存性はG20、端点との接続はE34・G25で明示的に検査する。

結果は [scalar_certificate.json](results/scalar_certificate.json)、
[localization_certificate.json](results/localization_certificate.json)、
[profile_certificate.json](results/profile_certificate.json)、
[resonant_gap_certificate.json](results/resonant_gap_certificate.json) に記録する。
これらは解析的推論やライブラリ・コード全体の独立した形式検証ではない。

### 11.2. 整数検算

[verify_exact.py](verification/verify_exact.py) は次を再実行した。

- (3.3) を $t=1,2,3,4,9$、$Q$ の次数350まで照合。
- Rogers–Ramanujan の和表示と積表示、および (2.3) を別計算で照合。
  5分解の別計算は Euler の五角数分子と分割数分母を使用し、$q$ の次数1754まで確認。
- $n=1,\ldots,256$ の全56,252,416係数を任意精度整数で検査。
  符号違反は0。対称性、剰余類和、最初の尾部公式も一致。
- $n=2,\ldots,256$ で零集合 (1.2) と完全一致。
- $n\le10$ は別の積の実装と全係数照合。
  $n=2,8,32,64,128,256$ は二つの素数法で $q=2$ の直接積と照合。

結果は [exact_certificate.json](results/exact_certificate.json)。
各 $n$ の全係数を十進整数と改行で連結した SHA-256 も記録した。
ハッシュは再実行時の一致確認用であり、それ自体が符号証明ではない。

### 11.3. 通常の多倍長診断

[diagnose_analytic.py](verification/diagnose_analytic.py) は80桁計算で
旧版の局所展開の48ケース、eta の正確な変換式 (5.3)–(5.4) の16ケースを照合した。
これらは有限個のサンプルであり、全角度・全パラメータの証明ではない。

さらに通常の分割数 $p(j)$ と累積 $P_c(K)=\sum_{j\le K}p(j)$ を使う

```math
d_K=-\sum_{b\in\mathbb Z}(-1)^b
P_c\left(K-\left\lfloor\frac{b(3b-1)}{10}\right\rfloor\right)
\tag{11.1}
```

を $K\le150$ で RR 級数による計算と照合した。
この公式は Euler の五角数指数が5を法として0、1、2のみになることと
$(gh_R+h_R^2-g^2)/(1-Q)$ から従う。

分割数を100000まで整数漸化式で計算し、(3.6) から得た十個の正確な係数を
端点主項と比較した。下表は剰余類3。剰余類4も結果ファイルに含む。

| $n$ | $m$ | $c_n(m)/[-P(n,m-5n)]$ |
|---:|---:|---:|
| 120000 | 900003 | 0.999590037443789303… |
| 160000 | 1200003 | 0.999645143026078… |
| 200000 | 1500003 | 0.999682716882507… |

全十例で (4.1) の条件を満たし、比は (4.3) の誤差範囲内にある。
巨大な $n$ の全係数を展開したものではない。
結果は [analytic_diagnostics.json](results/analytic_diagnostics.json)。

v0.3–v0.4 では [diagnose_threshold.py](verification/diagnose_threshold.py) により、
(8.12) の等比級数と直接和の60例、通常および重みを保持した局所剰余192例、
正規化した実際の積分8例も80桁で診断した。
積分診断は $\tau=0,8$、$n=10^7,2\cdot10^7$、根 $\ell=1,2$ を使い、
$u=\sqrt{nV}t\in[-8,8]$ の範囲で評価する。
全ての積分区間が $|t|\le0.1$ に含まれることもコードで確認する。
結果は [threshold_diagnostics.json](results/threshold_diagnostics.json)。
さらに従来の局所展開48例のうち小箱の32例に、積分型の剰余 (9.16) を適用し、
より小さい上界とも整合することを確認した。
v0.7 の [diagnose_grouped.py](verification/diagnose_grouped.py) は、有限積とB2補正を36例、
四根合成後の局所誤差・一次補正を120例、弱い剰余類の正規化積分を4例、80桁で照合する。
積分は $n=31147$、$\tau=0,11/2$、剰余類3・4を使い、一次EM補正も含めて比較する。
結果は [grouped_diagnostics.json](results/grouped_diagnostics.json)。
追加の [diagnose_localization.py](verification/diagnose_localization.py) は非共鳴部分集合の直接和76例、
等比級数核26例を80桁で照合する。結果は [localization_diagnostics.json](results/localization_diagnostics.json)。
端点の個別係数にはn=40000、K=39200の剰余類3・4を追加し、v0.6のw上限の範囲でも確認する。さらにn=31148、K=31147の二つの係数を追加し、w<=0.0013の範囲も照合した。
旧版の診断は、その旧版の適用範囲内で回帰検査として残した。
これらの有限診断を、§8–9 の全パラメータに対する解析的導出の代わりには使わない。

## 12. 限界、次の監査、閾値改善

本稿の証明案を検討するときは、次を優先する。

1. (2.5) を円周全体に適用し、(5.9) と eta 剰余の両方で $x_0$ を保持しているか。
2. (5.1) と (8.4) が全角度を覆い、局所診断だけに依存していないか。
3. 尾部展開が $x_0/w_0\to0$ を要求せず、共通指数を保存しているか。
4. 根の乗数、対数の枝、積分前因子、四根を合成する順番に誤りがないか。
5. 端のシフト付き鞍点と内部の未シフト鞍点の領域が §10 で隙間なく接続しているか。

これらの全推論を Lean 等で形式検証したわけではない。
外部の専門家への連絡や査読を完了条件にはしていないが、
その分、公開時にも独立監査前という状態を明示する。

v0.2 の $1/\sqrt n$ 型の誤差は v0.3 で $1/n$ 型に改善した。
v0.4 は同じ係数位置で分散・高次モーメント・位相を比較し、
積分型の有限和剰余と全角度の定数改善を加えた。
これにより v0.3 の一律比較が要求した約 $1.6961\cdot10^9$ という制約も取り除いた。

v0.7は複数のパラメータを同時探索した後、$N=31147$ で区間証明書を再生成した。
端点誤差の上限は約0.84662267で、採用した0.85未満に収まる。
主円弧の誤差・位相比の上界も1未満である。

今回「限界」と確認した対象は、固定した $w_{\max}=0.0013$ における端点条件である。

```math
31146<\frac{C}{25w_{\max}^2}=31146.6805557\ldots<31147.
\tag{12.1}
```

[threshold_boundary_certificate.json](results/threshold_boundary_certificate.json) は両側を区間検証する。
従ってこの端点条件のまま整数閾値を31146以下へ下げることはできない。
これは全てのパラメータや証明の組み替えに対する最適性の証明ではない。
$w_{\max}$ をさらに変える、端点の誤差式を改善する、別の係数領域の被覆を使うなどの変更は未排除。
探索の最良値と数学的な絶対下限を区別する。

閾値の履歴は $10^{16}\to4\cdot10^{13}\to2\cdot10^9\to5\cdot10^6\to300000\to40000\to31147$。
$257\le n<31147$ の全体、全$n$の零集合分類、解析的証明の独立監査は引き続き未了。
この残りを全て総当たりできる実用的な方法を得たわけではない。

## 13. AI 使用、来歴、参考文献

本稿は AI との対話から作成された研究作業ノートを基に、
OpenAI Codex が構成を統合し、説明を補い、検証コードを新規作成して実行したものである。
v0.2–v0.7 の閾値改良に用いた追加の補題、定数評価、区間証明書も Codex が導出・作成した。
v0.5–v0.7 はユーザーが提供した別のAIによる改善提案・監査意見を出発点とし、Codex が再導出・実装・区間検証・検算を行った。
v0.6では平滑化に関する反論も照合し、半径を分子・分母の両方に保持する導出を本文へ補足した。
共同ソートの全分母列挙や端点の奇数次消去は今回の証明には使用していない。
人間による全証明の独立検証を受けたとは主張しない。
原ノートに記録された過去の実行履歴は今回確認しておらず、
ここに添付する数値結果は今回の再計算による。
著者名は記載しない。趣味の研究記録として保存・公開する。

原ノート自体は本リポジトリには収録していない。
本稿と検証コード・出力の SHA-256 一覧を添付する。
本稿は先行する作業ノートを定理の外部根拠として引用せず、
必要な推論を本文に再記した。

- **[RR]** NIST DLMF, [§17.2(vi), Rogers–Ramanujan identities](https://dlmf.nist.gov/17.2#vi).
- **[Eta]** NIST DLMF, [§23.18, Dedekind eta modular transformations](https://dlmf.nist.gov/23.18).
- **[Coth]** NIST DLMF, [§4.36.3, partial fraction expansion of coth](https://dlmf.nist.gov/4.36#E3).
- **[W21]** Liuquan Wang, [Sign Changes of Coefficients of Powers of the Infinite Borwein Product](https://arxiv.org/html/2108.03932v3), Proposition 3.6 とその証明.
- **[WK22]** Chen Wang and Christian Krattenthaler, [An asymptotic approach to Borwein-type sign pattern theorems](https://arxiv.org/html/2201.12415v1), 特に §11 の第三予想の相殺に関する議論.

RR 恒等式、無限積の5分解、eta 変換、鞍点法は既知の道具である。
非負核恒等式、有限尾部への適用、eventual theorem、有効閾値、零集合分類の
文献上の新規性については、網羅的調査を完了していない。
