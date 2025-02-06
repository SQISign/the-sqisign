Signature size and compression methods
======================================

We describe the signature compression, and discuss some tradeoffs.
Concrete numbers are given for NIST level 1, λ=128b=16B

## Basic signature

Parameters: p ~ 2 λ bits; p=251b=32B

The signature consists of:
- j(E_aux): 2 log_2(p) ~ 4 λ bits; j=502b=64B
- the challenge scalar c: λ bits; c=128b=16B
- a representation of P_chl, Q_chl given by a matrix M in terms of a canonical basis of E_aux. These points are of order dividing 2^{e_rsp}, with e_rsp ~ λ.
  By Minkowski, we are guaranteed a response of size $\le 2^{3/2}/π p^{1/2} ~ 0.9 p^{1/2}$ ~ λ bits. The matrix thus needs 4λ bits.
  M=512b=64B

Total: 8λ bits.

## Hints

We currently add some hints for a faster verification:
- 16b = 2B for the canonical basis generation on E_aux
  Having it speed up the canonical basis step.

- n_bt: 2^{n_bt} is the amount of backtracking the response isogeny does through the commitment.
  Having it allows to not compute useless isogenies in the challenge that
  will get backtracked through the response; so save a few dimension 1
  2-isogenies computations.

- r': 2-adic valuation of the degree of the response.
  The (non backtracking part of the) response is of degree q=2^{r'} q'
  The part 2^{r'} is computed through a dim 1 isogeny.
  The part q' is computed through a dim 2 isogeny.
  We could recover r' directly by computing the order of the Weil pairing of P_chl, Q_chl (formally if they are given by matrix coefficients)

-> n_bt, r': log(λ) bits; 2x7b = 14b=2B

Total with hints: 9λ+2+(ε=2log(λ)); 1172b=146.5B

Remark 1: for the heuristic version, n_bt is not needed because we sample the response to be non backtracking.

Remark 2: for the rigorous version, the first thing we compute is E'_aux and
a representation of the isogeny E_com x E'_aux -> E'_chl x E_aux.
But for the signature, we would need to output j(E_com), j(E'_aux), adding
4λ to the signature size. That's why we compute this isogeny explicitly in the signature, to get a representation of the dual isogeny
  E'_chl x E_aux -> E_com x E'_aux.
Here E'_chl can be recovered from the challenge data, so we only need
to put j(E_aux) in the signature. So we compress by 4λ compared to the
representation above, at the cost of an extra dim 2 isogeny.

## Pairings

We could save λ by only giving 3 out of 4 of the coefficients of the matrix
M, which would then use 3λ bits rather than 4λ bits. The last coefficient
can be recovered by a pairing, using the fact that the kernel of the dim 2
isogeny is isotropic.

In that version, we really need to output n_bt and r' (or at least n_bt+r').

## Misc optimisation

A dim 1 isogeny of degree 2^r only needs r bits to be represented (r+1
actually), a dim 2 isogeny of type (2^r, 2^r) needs 3r bits.
This means that the higher n_bt and r' are, the less space we need for the
matrix M; so we could use a variable bit length encoding to represent them.

In practice, this would work as follow: the matrix M encodes P_chl, Q_chl
of different order. P_chl is encoded with respect to 2^{n_bt} times the
canonical basis, while Q is encoded with respect to 2^{n_bt+r′} times the
canonical basis (and, if needed, only by one coefficient rather than two,
the second one can be found by a pairing.)

For P, the r'-most significant bits of its coefficients a,b encode the
kernel of the dim 1 isogeny, while the remaining bits encode what we need
for the dim 2 representation. We can scale these most significant bits to
encode the dim 1 isogeny, so that they are of the form (1+2^{r'} ...)
P_canonical, (c+2^{r'} ...) Q_canonial.
(+ a flag in case the kernel is given by 2^{e-r'} Q_canonical)

This allows to encode n_{bt}, r' in such a way to add 2b rather than 2B to
the signature.

## Max compression

Using pairings to save λ, removing the 2B hints for the canonical basis
generation, and using the trick above to gain 2B-2b for n_{bt} and r' (or
even removing them entirely), the signature can be compressed to 8λ.

## Size/speed tradeoff

In the other direction, one could give the $x$-coordinates of P_chl, Q_chl
instead of the matrix M. This would take 2*4λ=8λ, while M only takes 4λ or
3λ with pairings. But we would skip the matrix action evaluation, which can
take as much as 25% of the verification time.
