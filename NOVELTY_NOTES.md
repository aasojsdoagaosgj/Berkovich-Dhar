# Novelty / prior-art notes

**Search date:** 18 September 2026.

This file records a targeted literature check. It is not a formal priority determination and cannot guarantee that no unpublished, unindexed, or differently worded result exists.

## What was already in the literature

Berkovich and Dhar's *New Borwein-Type Conjectures* (arXiv:2407.13788; *Experimental Mathematics*) explicitly formulates five modulus-3 cases with powers $4,5,6,7,8$ and two modulus-5 cases with powers $2,3$. The present repository concerns the six cases with powers $5,6,7,8$ modulo 3 and $2,3$ modulo 5; the quartic modulus-3 case is retained only as an analytic dependency.

Wang and Krattenthaler's arXiv:2201.12415 develops the asymptotic/saddle framework for Borwein-type sign theorems and is the principal methodological precursor. Liuquan Wang's arXiv:2108.03932 supplies asymptotic sign information for powers of the corresponding infinite products and is used here for the initial-degree input.

Thus the broad idea "use asymptotic/saddle methods" is **not** new.

## What appears to be new in this repository

A targeted search did not locate a prior proof of the following package of results:

1. Eventual sign-pattern theorems for all six cases
   $(3,5),(3,6),(3,7),(3,8),(5,2),(5,3)$.
2. Exact monotone phase equations determining the transition limits in these cases.
3. The refinement
```math
   k_n=\alpha n^2+\beta n+O(1)
```
   for each transition.
4. The modulus-5 four-root factorisation used to prove uniqueness of both relevant transitions.
5. The dedicated treatment of the modulus-5 cube's vanishing infinite-product residues, including the coefficientwise first-tail kernel and the whole-circle subtraction preserving the small tail factor.
6. Correction of several decimal transition limits printed in the conjecture paper; the repository's values come from the exact phase equations rather than extrapolation from finite $n$.

## Searches checked

The search included the exact paper/conjecture names, arXiv identifiers, author publication listings, combinations of "Berkovich Dhar" with "Conjecture 2.1", "Conjecture 2.2", "Borwein sign pattern", "asymptotic", and the relevant modulus/power pairs, as well as searches for the high-precision transition constants.

At the search date, the journal landing page for *New Borwein-Type Conjectures* showed no CrossRef citations, and the authors' indexed/publication listings did not surface a later resolution of these six cases. These are useful indicators but not proofs of novelty.

## Conservative public wording

A safe description is:

> This repository gives a proof draft of the eventual (large-$n$) form of six Berkovich–Dhar Borwein-type cases, using the Wang–Krattenthaler asymptotic framework together with additional transition and endpoint analysis. A targeted literature search on 18 September 2026 did not locate a prior proof of these six eventual theorems.

Avoid describing this repository as a proof of the full Berkovich–Dhar conjectures, since the original conjectures quantify over every positive integer $n$.
