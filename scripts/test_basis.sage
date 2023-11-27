#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')

################################################################

from parameters import p, B, f, Tpls, Tmin, Dcom, Dchall
T = Tpls * Tmin

################################################################

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')

Fp2.<i> = GF((p,2), modulus=[1,0,1])
Fp4 = Fp2.extension(2,'u')
E = EllipticCurve(Fp4, [1,0])
assert E.j_invariant() == 1728
assert E.is_supersingular()
assert E.change_ring(Fp2).frobenius() == -p
assert E.order() == (p^2-1)^2

from sage.groups.generic import order_from_multiple
x = Fp4.gen()
while True:
    x += 1
    try:
        P = E.lift_x(x)
    except ValueError:
        continue
    o = order_from_multiple(P, p^2-1)
    if (T<<f).divides(o):
        P *= o // (T<<f)
        P.set_order(T<<f)
        break

x = Fp4.gen()
while True:
    x += 1
    try:
        Q = E.lift_x(x)
    except ValueError:
        continue
    o = order_from_multiple(Q, p^2-1)
    if not (T<<f).divides(o):
        continue
    Q *= o // (T<<f)
    Q.set_order(T<<f)
    if order_from_multiple(P.weil_pairing(Q, T<<f), T<<f, operation='*') == T<<f:
        break


bases = {
        'A': Tpls,
        'B': Tmin,
        '2': 2^f,
        '3': 3^valuation(Tpls,3),
    }

from cformat import EcPointX, Object, ObjectFormatter
objs = ObjectFormatter(
    [Object('const fp2_t', f'xP{k}', EcPointX(ZZ(P.order()/v)*P, Fp2))
    for k,v in bases.items()] +
    [Object('const fp2_t', f'xQ{k}', EcPointX(ZZ(Q.order()/v)*Q, Fp2))
    for k,v in bases.items()] +
    [Object('const fp2_t', f'xPQ{k}', EcPointX(ZZ((P-Q).order()/v)*(P-Q), Fp2))
    for k,v in bases.items()]
    )

with open('test-basis.h', 'w') as hfile:
    print('#ifndef TEST_BASIS_H', file=hfile)
    print('#define TEST_BASIS_H', file=hfile)
    print('#include \"fp2.h\"', file=hfile)
    objs.implementation(file=hfile)
    print('#endif', file=hfile)