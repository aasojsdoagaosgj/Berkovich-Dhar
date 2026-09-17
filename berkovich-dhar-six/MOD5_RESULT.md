# The two eventual modulus-5 theorems

18 September 2026. Anonymous hobby research manuscript. Mathematical
derivation, code and checks by OpenAI Codex. Not independently reviewed
or formally verified; no publication or external contact has occurred.

For s=2,3, define c_j(n)=[q^j]P_{n,5}(q)^s. The degree is 10*s*n^2.
There is an integer N_s such that for all n>=N_s the following holds:

* c_{5k}(n)>0 for every index in its degree range;
* in each of residues 3 and 4 there is an index k_{n,s,a} such that
  c_{5k+a}(n)>=0 before and at that index and <=0 afterwards;
* k_{n,s,a}/n^2 converges to alpha_{s,a} defined by (P2),(P7) in
  MOD5_COMPACT.md, and in fact
  k_{n,s,a}=alpha_{s,a}*n^2+beta_{s,a}*n+O(1).

Choose k_{n,s,a} as one less than the first index with negative
coefficient. This definition handles the possible zero at the
transition and the terminal zero block in the cube's residue 3.
It is NOT the last nonnegative coefficient of the entire polynomial.

The exact constants are -s*R'(tau), where the unique tau>0 solves
V(e^{-tau})=-pi/(10*s) for a=3 and U(e^{-tau})=-pi/(10*s) for a=4.
Numerical evaluations (not rational interval certificates) are

| s | residue 3 limit | residue 4 limit |
|---|---|---|
| 2 | 1.1354609594225800624... | 1.0052808235905823023... |
| 3 | 2.1368999196948714797... | 1.9842962321701810724... |

These differ from the decimal estimates in
[Berkovich--Dhar, Conjecture 2.2](https://arxiv.org/html/2407.13788v3#S2).
The proof concerns its qualitative sign assertions and convergence,
with the limits defined exactly above. No numerical threshold N_s
or all-positive-n theorem is asserted.

## Proof of coverage and reciprocity

Initial degrees through 5n agree with the infinite products, whose
signs are proved in INFINITE_SIGNS.md. For s=2 they are +,-,-,+,+,
with just the exceptional zero of degree 9. For s=3 they are
+,-,0,+,0, with residues 2 and 4 identically zero.

For s=3, MOD5_FIRST_TAIL.md proves the strict signs of those two
weak residues from degrees 5n+2 and 5n+4, respectively, through
the entire first tail below 10n. This is an all-n algebraic identity
and a coefficientwise positive-kernel proof.

MOD5_ENDPOINT.md covers all remaining lower-half degrees up to an
interval that overlaps the compact range. Its weak case uses the
shifted degree j-5n and the exact subtraction of G_5^3 before every
error estimate; its large parameter N exceeds sqrt(n)/2. In this
endpoint region residues 0,3,4 are positive and residues 1,2 negative.

MOD5_COMPACT.md covers 5n^2*(-H'(11/2))<=j<=5sn^2. It shows
that residues 0,1,2 retain these signs and that residues 3 and 4 each
change from positive to negative exactly once. At most one zero
can occur near either transition, by the normalized adjacent-degree
comparison. The endpoint upper boundary has -H'(5)>-H'(11/2),
so every integer degree in the lower half is covered.

There are 4*n*s reciprocal factors with total degree 10*s*n^2.
Thus c_j=c_{10*s*n^2-j}. Residue 0 reflects to itself, 3 reflects
to 2 and 4 to 1. In the lower half residues 1 and 2 are nonpositive
everywhere. Therefore the upper-half residues 3 and 4 are
nonpositive everywhere and cannot produce a second positive block.
This proves both whole-polynomial sign assertions.

For the square, degree 9 is a persistent zero in the initial positive
block of residue 4, so the theorem does not claim a unique zero
globally. For the cube, the initial zeros of residue 2 reflect to
terminal zeros of residue 3. Initial zeros of residue 4 also remain.
These are compatible with the conjecture's weak inequalities.

Finally the compact expansion yields (P7), hence both limits and
linear corrections. All sufficiently-large conditions are uniform,
and only finitely many are used, establishing a common threshold
for the two powers. This proves the stated eventual theorems,
subject to the written analytic dependencies, not merely to the
finite numerical samples.
