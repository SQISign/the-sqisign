#!/usr/bin/env sage
proof.all(False)  # faster

from sage.misc.banner import require_version
if not require_version(10, 0, print_message=True):
    exit('')

################################################################

from parameters import p, f
from torsion_basis import even_torsion_basis_E0

################################################################

from sage.groups.generic import order_from_multiple
pari.allocatemem(1 << 34)  # 16G

if p % 4 != 3:
    raise NotImplementedError('requires p ≡ 3 (mod 4)')

assert (1 << f).divides(p + 1)
Fp2.<i> = GF((p,2), modulus=[1,0,1])

sqrtm1 = min(Fp2(-1).sqrt(all=True))

def compute(q, mat, idl, iso1q):
    print(f'\x1b[33m{q = }\x1b[0m')
    E0 = EllipticCurve(Fp2, [1,0])
    E0.set_order((p+1)^2)

    if q == 1:
        E1 = E0
        P1, Q1 = even_torsion_basis_E0(E1, f)
        print(f'E0 = {E1}')
        print(f'P0 = {P1}')
        print(f'Q0 = {Q1}')

    else:
        Quat.<i,j,k> = QuaternionAlgebra(-1, -p)
        I = Quat.ideal(map(Quat, idl))
#        print(f'{I = }')
        O0 = Quat.quaternion_order(list(map(Quat, orders[0][2])))
#        print(f'{O0 = }')
        O1 = I.right_order()
#        print(f'{O1 = }')
        assert I.left_order() == O0
        assert O0.is_maximal() and O1.is_maximal()
        assert I.norm() % 2

        from deuring2d import Deuring2D
        ctx = Deuring2D(p)
        assert ctx.O0.order == O0
        assert ctx.E0 == E0
        ctx.sqrtm1 = sqrtm1

        P0, Q0 = data[0][1]

        for deg in range(1,10):

            print(f'trying {deg = }...')
            ctx.e = E0.cardinality(extension_degree=2^deg).sqrt().valuation(2) - 1

            first = True
            for suitable in ctx.SuitableIdeals(I, attempts=10**6, bound=10**3):

                if first:
                    Fbig.<U> = Fp2.extension(2^deg)
                    ctx.E0 = E0.change_ring(Fbig)
                    ctx.P = P0.change_ring(Fbig)
                    ctx.Q = Q0.change_ring(Fbig)
                    assert ctx.e == ctx.E0.order().sqrt().valuation(2) - 1
                    for _ in range(ctx.e - f):
                        ctx.P = ctx.P.division_points(2)[0]
                        ctx.Q = ctx.Q.division_points(2)[0]
                    ctx.P.set_order(multiple=2^ctx.e)
                    ctx.Q.set_order(multiple=2^ctx.e)
                first = False

                try:
                    E1, P1, Q1 = ctx.IdealToIsogeny(I, suitable=suitable)
                    break
                except Deuring2D.Failure:
                    continue

            else:
                continue
            break

        else:
            raise NotImplementedError('Deuring2D failed')

        E1 = E1.change_ring(Fp2)

        j = GF(p)(E1.j_invariant())
        X = polygen(GF(p))
        for A,_ in sorted((256*(X^2-3)^3 - (X^2-4)*j).roots()):
            E1_ = EllipticCurve(Fp2, [0,A,0,1,0])
            try:
                iso = min(E1.isomorphisms(E1_))
                break
            except ValueError:
                pass
        E1 = iso.codomain()
        P1 = iso._eval(P1)
        Q1 = iso._eval(Q1)
        print(f'{E1 = }')

        P1 *= ctx.P.order() // P0.order()
        Q1 *= ctx.Q.order() // Q0.order()
        P1 = P1.change_ring(Fp2)
        Q1 = Q1.change_ring(Fp2)
        print(f'{P1 = }')
        print(f'{Q1 = }')
        P1.set_order(P0.order())
        Q1.set_order(Q0.order())
        assert P0.order() == Q0.order() == P1.order() == Q1.order() == 2^f

        assert P1.weil_pairing(Q1,2^f) == P0.weil_pairing(Q0,2^f)^I.norm()

    if q == 1:
        endo_i, = (a for a in E1.automorphisms() if a.scaling_factor() == sqrtm1)
    else:
        iso = E1.isomorphism(min(Fp2(-q).sqrt(all=True)), is_codomain=True)
        try:
            endo_i = iso * E1.isogeny(None, codomain=iso.domain(), degree=q)
        except ValueError:
            assert False
#    assert endo_i^2 == -q

    endo_1 = E1.scalar_multiplication(1)
    endo_j = E1.frobenius_isogeny()
    endo_k = endo_i * endo_j

    if __debug__:
        R = E1.random_point()
        assert (endo_i^2)(R) == -q*R
        assert (endo_j^2)(R) == -p*R
        assert (endo_j*endo_i)(R) == -(endo_i*endo_j)(R)

    denom = mat.denominator()
    coprime = denom.prime_to_m_part(lcm(P1.order(), Q1.order()))
    P1d, Q1d = (inverse_mod(coprime, T.order()) * T for T in (P1, Q1))

    denom //= coprime

    extdeg = next(d for d in range(1,denom+1) if ((denom<<f)^2).divides(E1.order(extension_degree=d)))
    if extdeg == 1:
        Fbig = Fp2
    else:
        Fbig.<U> = Fp2.extension(extdeg)

    P1d, Q1d = (T.change_ring(Fbig) for T in (P1d, Q1d))
    P1d.set_order(multiple=denom<<f)
    for l,m in denom.factor():
        for i in range(m):
            assert l.divides(P1d.order())
            P1d = P1d.division_points(l)[0]
            P1d.set_order(multiple=denom<<f)
            for Q1d_ in Q1d.division_points(l):
                o = order_from_multiple(P1d.weil_pairing(Q1d_, P1d.order()), denom<<f, operation='*')
                if o == P1d.order():
                    Q1d = Q1d_
                    break
            else:
                assert False
    assert hasattr(P1d, '_order')
    Q1d.set_order(multiple=denom<<f)

    denom *= coprime

    PQ1d = P1d, Q1d
#    mat1 = matrix(Zmod(1<<f), [endo_1._eval(T).log(PQ1d) for T in PQ1d])
#    assert mat1 == 1            # identity; omit
    mati = matrix(Zmod(1<<f), [endo_i._eval(T).log(PQ1d) for T in PQ1d])
    matj = matrix(Zmod(1<<f), [endo_j._eval(T).log(PQ1d) for T in PQ1d])
#    matk = matrix(Zmod(1<<f), [endo_k._eval(T).log(PQ1d) for T in PQ1d])
#    assert matk == matj * mati  # redundant; omit
    matk = matj * mati

    gens = []
    for row in denom * mat:
        endo = sum(ZZ(c)*e for c,e in zip(row, (endo_1,endo_i,endo_j,endo_k)))
        gens.append(endo)
    gen1, gen2, gen3, gen4 = gens

    assert mat[0] == vector((1,0,0,0))
#    mat1 = matrix(ZZ, [gen1._eval(T).log(PQ1d) for T in PQ1d]) / denom
#    assert mat1 == 1            # identity; omit
    mat2 = matrix(ZZ, [gen2._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat3 = matrix(ZZ, [gen3._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat4 = matrix(ZZ, [gen4._eval(T).log(PQ1d) for T in PQ1d]) / denom
    mat2, mat3, mat4 = (M.change_ring(Zmod(1<<f)) for M in (mat2,mat3,mat4))

    A = E1.a2()
    assert E1.a_invariants() == (0,A,0,1,0)

    return (A, (A+2)/4), (P1, Q1), (mati,matj,matk), (mat2,mat3,mat4)

################################################################

from maxorders import orders

print('qs:', [q for q,_,_,_,_,_ in orders])

todo = [(q, mat*iso1q, idl, iso1q) for q,iso1q,mat,_,idl,_ in orders]
data = [None] * len(todo)

assert todo[0][0] == 1
data[0] = compute(*todo[0])  # compute this first; we need it for the others
print(f'[\x1b[32m+\x1b[0m] finished precomputation for \x1b[36mq = {todo[0][0]}\x1b[0m.')

####XXX
##for idx,inp in enumerate(todo[1:],1):
##    data[idx] = compute(*inp)
##    print(f'[\x1b[32m+\x1b[0m] finished precomputation for \x1b[36mq = {inp[0]}\x1b[0m.')
##todo = []
####XXX

for (inp,_),res in parallel(8)(compute)(todo[1:]):
    q,_,_,_ = inp
    idx, = (i for i,(qq,_,_,_) in enumerate(todo) if qq == q)
    assert data[idx] is None
    data[idx] = res
    print(f'[\x1b[32m+\x1b[0m] finished precomputation for \x1b[36m{q = }\x1b[0m.')

################################################################

from cformat import FpEl, Ibz, Object, ObjectFormatter

def Fp2_to_list(el):
    return [FpEl(int(c), p, True) for c in Fp2(el)]

def basis2field(P, Q):
    vs = [
            [Fp2_to_list(T[0]), Fp2_to_list(T[2])]
            for T in (P,Q,P-Q)
        ]
    return vs

################################################################

objs = ObjectFormatter([
        Object('curve_with_endomorphism_ring_t[]', 'CURVES_WITH_ENDOMORPHISMS',
            [
                [
                    [Fp2_to_list(A), Fp2_to_list(1),                    # ec_curve_t A, C
                     [Fp2_to_list(A24), Fp2_to_list(1)], "true"],       # ec_curve_t A24, is_A24_computed_and_normalized
                    basis2field(*basis),                                # ec_basis_t
                    [[Ibz(v) for v in vs] for vs in mati.transpose()],  # ibz_mat_2x2_t
                    [[Ibz(v) for v in vs] for vs in matj.transpose()],  # ibz_mat_2x2_t
                    [[Ibz(v) for v in vs] for vs in matk.transpose()],  # ibz_mat_2x2_t
                    [[Ibz(v) for v in vs] for vs in mat2.transpose()],  # ibz_mat_2x2_t
                    [[Ibz(v) for v in vs] for vs in mat3.transpose()],  # ibz_mat_2x2_t
                    [[Ibz(v) for v in vs] for vs in mat4.transpose()],  # ibz_mat_2x2_t
                ]
                for (A,A24),basis,(mati,matj,matk),(mat2,mat3,mat4)
                in data
            ])
    ])

with open('include/endomorphism_action.h','w') as hfile:
    with open('endomorphism_action.c','w') as cfile:
        print(f'#ifndef ENDOMORPHISM_ACTION_H', file=hfile)
        print(f'#define ENDOMORPHISM_ACTION_H', file=hfile)
        print(f'#include <sqisign_namespace.h>', file=hfile)
        print(f'#include <ec.h>', file=hfile)
        print(f'#include <quaternion.h>', file=hfile)
        print(f'#include <stddef.h>', file=cfile)
        print(f'#include <stdint.h>', file=cfile)
        print(f'#include <endomorphism_action.h>', file=cfile)

        print('''
/** Type for precomputed endomorphism rings applied to precomputed torsion bases.
 *
 * Precomputed by the precompute scripts.
 *
 * @typedef curve_with_endomorphism_ring_t
 *
 * @struct curve_with_endomorphism_ring
 **/
typedef struct curve_with_endomorphism_ring {
    ec_curve_t curve;
    ec_basis_t basis_even;
    ibz_mat_2x2_t action_i, action_j, action_k;
    ibz_mat_2x2_t action_gen2, action_gen3, action_gen4;
} curve_with_endomorphism_ring_t;
              '''.strip(), file=hfile)

        print(f'#define CURVE_E0 (CURVES_WITH_ENDOMORPHISMS->curve)', file=hfile)
        print(f'#define BASIS_EVEN (CURVES_WITH_ENDOMORPHISMS->basis_even)', file=hfile)
        print(f'#define ACTION_I (CURVES_WITH_ENDOMORPHISMS->action_i)', file=hfile)
        print(f'#define ACTION_J (CURVES_WITH_ENDOMORPHISMS->action_j)', file=hfile)
        print(f'#define ACTION_K (CURVES_WITH_ENDOMORPHISMS->action_k)', file=hfile)
        print(f'#define ACTION_GEN2 (CURVES_WITH_ENDOMORPHISMS->action_gen2)', file=hfile)
        print(f'#define ACTION_GEN3 (CURVES_WITH_ENDOMORPHISMS->action_gen3)', file=hfile)
        print(f'#define ACTION_GEN4 (CURVES_WITH_ENDOMORPHISMS->action_gen4)', file=hfile)
        print(f'#define NUM_ALTERNATE_STARTING_CURVES {len(data)-1}', file=hfile)
        print(f'#define ALTERNATE_STARTING_CURVES (CURVES_WITH_ENDOMORPHISMS+1)', file=hfile)

        objs.header(file=hfile)
        objs.implementation(file=cfile)

        print(f'#endif', file=hfile)
