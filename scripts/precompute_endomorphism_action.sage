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

endo_1 = E.scalar_multiplication(1)
endo_i = E.automorphisms()[-1]
endo_j = E.frobenius_isogeny()
endo_k = endo_i * endo_j

if 0:  # skipped for speed, for now
    assert endo_i^2 == E.scalar_multiplication(-1)
    assert endo_j^2 == E.scalar_multiplication(-p)
    assert endo_j * endo_i == - endo_i * endo_j
else:
    R = E.random_point()
    assert (endo_i^2)(R) == -1*R
    assert (endo_j^2)(R) == -p*R
    assert (endo_j*endo_i)(R) == -(endo_i*endo_j)(R)

def half_endo(summands):
    def _eval(P):
        E = P.curve()
        assert P in E
        F = E.base_field()
        if (halves := P.division_points(2)):
            Q = halves[0]
        else:
            Q = E.change_ring(F.extension(2,'v'))(P)
        R = sum(endo._eval(Q) for endo in summands)
        return E(R)
    return _eval

gen1 = endo_1._eval
gen2 = endo_i._eval
gen3 = half_endo([endo_i, endo_j])
gen4 = half_endo([endo_1, endo_k])

################################################################

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

def dlp(P, Q, R):
    n = P.order()
    assert P.order() == Q.order()
    assert R.order().divides(P.order())
    e = Fp2(P.weil_pairing(Q, n))
    a = Fp2(R.weil_pairing(Q, n)).log(e)
    b = Fp2(P.weil_pairing(R, n)).log(e)
    assert a*P + b*Q == R
    return a, b

def matrix_of_isogeny(phi):
    imP, imQ = map(phi, (P,Q))
    vecP = dlp(P, Q, imP)
    vecQ = dlp(P, Q, imQ)
    mat = matrix(Zmod(T<<f), [vecP, vecQ]).transpose()
    assert imP == ZZ(mat[0][0])*P + ZZ(mat[1][0])*Q
    assert imQ == ZZ(mat[0][1])*P + ZZ(mat[1][1])*Q
    return mat

#mat1 = matrix_of_isogeny(endo_1)
mati = matrix_of_isogeny(endo_i)
matj = matrix_of_isogeny(endo_j)
matk = matrix_of_isogeny(endo_k)
#assert mat1 == 1    # identity; omit

#mat1 = matrix_of_isogeny(gen1)
mat2 = matrix_of_isogeny(gen2)
mat3 = matrix_of_isogeny(gen3)
mat4 = matrix_of_isogeny(gen4)
#assert mat1 == 1    # identity; omit

################################################################

Quat.<i,j,k> = QuaternionAlgebra(-1, -p)
O0 = Quat.quaternion_order([1, i, (i+j)/2, (1+k)/2])

assert Dcom % 2 == 1  # odd
mat = block_matrix(Zmod(Dcom), [[identity_matrix(2), mati, matj, matk]])[:,::2]
ker = list(map(Quat, mat.right_kernel_matrix()))
idealP = sum((O0*g for g in ker), O0*Dcom)
assert idealP.norm() == Dcom
for b in idealP.basis():
    assert sum(Mod(c,Dcom)*g for c,g in zip(b,(1,mati,matj,matk)))[:,0] == 0  # kills P
for v in (ZZ^4):
    idealPgen = sum(c*g for c,g in zip(v, idealP.basis()))
    if vector(list(idealPgen)).denominator() == 2:
        idealPgen *= 2
    if gcd(idealPgen.reduced_norm(), Dcom^2) == Dcom:
        break
assert idealP == O0*Dcom + O0*idealPgen

mat = mat   # still
rhs = vector(Zmod(Dcom), [0,1])
cs = mat.solve_right(rhs)
distorter = Quat(cs)
assert sum(Mod(c,Dcom)*g for c,g in zip(distorter,(1,mati,matj,matk))).columns()[0] == vector((0,1))  # maps P->Q

################################################################

from cformat import Ibz, Basis, Object, ObjectFormatter

# def field2limbs(el):
#     l = 1 + floor(log(p, 2**64))
#     el = Fp2(el)
#     vs = [[(int(c) >> 64*i) % 2**64 for i in range(l)] for c in el]
#     return vs

# def fmt_basis(name, P, Q):
#     vs = [
#             [field2limbs(T[0]), field2limbs(T[2])]
#             for T in (P,Q,P-Q)
#         ]
#     return Object('ec_basis_t', name, vs)

bases = {
        'EVEN': 1<<f,
        'ODD_PLUS': Tpls,
        'ODD_MINUS': Tmin,
        'COMMITMENT_PLUS': gcd(Tpls, Dcom),
        'COMMITMENT_MINUS': gcd(Tmin, Dcom),
        'CHALLENGE': Dchall,
    }

assert P.order() == Q.order()

objs = ObjectFormatter([
        Object('ec_basis_t', f'BASIS_{k}',
            Basis(p, Fp2, ZZ(P.order()/v)*P, ZZ(Q.order()/v)*Q)
        )
            for k,v in bases.items()
    ] + [
        Object('ec_curve_t', 'CURVE_E0', [[[int(0)]], [[int(1)]]]),
        Object('ec_point_t', 'CURVE_E0_A24', [[[int(0)]], [[int(1)]]]),
        Object('ibz_mat_2x2_t', 'ACTION_I', [[Ibz(v) for v in vs] for vs in mati]),
        Object('ibz_mat_2x2_t', 'ACTION_J', [[Ibz(v) for v in vs] for vs in matj]),
        Object('ibz_mat_2x2_t', 'ACTION_K', [[Ibz(v) for v in vs] for vs in matk]),
        Object('ibz_mat_2x2_t', 'ACTION_GEN2', [[Ibz(v) for v in vs] for vs in mat2]),
        Object('ibz_mat_2x2_t', 'ACTION_GEN3', [[Ibz(v) for v in vs] for vs in mat3]),
        Object('ibz_mat_2x2_t', 'ACTION_GEN4', [[Ibz(v) for v in vs] for vs in mat4]),
        Object('quat_alg_elem_t', 'COMMITMENT_IDEAL_UNDISTORTED_GEN', [Ibz(1), [Ibz(ZZ(v)) for v in idealPgen]]),
        Object('quat_alg_elem_t', 'COMMITMENT_IDEAL_DISTORTION_ENDO', [Ibz(1), [Ibz(ZZ(v)) for v in distorter]]),
    ])

with open('include/endomorphism_action.h','w') as hfile:
    with open('endomorphism_action.c','w') as cfile:
        print(f'#include <intbig.h>', file=hfile)
        print(f'#include <ec.h>', file=hfile)
        print(f'#include <quaternion.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <endomorphism_action.h>', file=cfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

