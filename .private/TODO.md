A list of short to long term possible optimisations.

# Short term

- Precompute the dim 1 and dim 2 isogeny strategies
- Precomputing the gluing isogenies from $E_0$.
  We compute several dim 2 isogenies from E_0 x E_0. The gluing step is
  quite costly, and we could simply precompute all 6 possible gluings from
  E_0.
- At several points (when our points have full 2^e order), we could use the Tate pairing instead of the Weil pairing.
- Pairings involve arithmetic doublings, which we could reuse later.
- For the different arithmetic operations on basis we need: investigate whether it is faster to use x(P), y(P), x(Q), y(Q) (in Edwards coordinates) or x(P), x(Q), x(P+Q).
- Faster low level arithmetic for F_{p^2}

# Medium term

Any improvement on dim 2 isogenies would improve the signature and verification:
- use theta squared coordinates rather than theta coordinates, this would save some square roots when we don't have the extra torsion available
- directly compute the gluing isogeny from Montgomery coordinates on the curve to theta on the surface without going through the conversion step of theta on the curve
- likewise: directly compute the splitting isogeny to Montgomery coordinates on the curves

For IdealToIsogeny: investigate if we can find a better way to sample quaternion elements to find one of u, v a sum of two squares and save one dim 2 isogeny.

# Long term

- Split 2^n dim 2 isogenies into 4-isogenies rather than 2-isogenies
- Investigate the use of Weil's restriction (cf Costello, Reijnders, Santos)
- Investigate the use of challenge by torsion rather than challenge by isogeny
