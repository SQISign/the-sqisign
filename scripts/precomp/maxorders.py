from sage.all import *

from sage.misc.banner import require_version
if not require_version(10, 5, print_message=True):
    exit('')

from parameters import p, num_orders as num

################################################################

# Underlying theory:
# - Ibukiyama, On maximal orders of division quaternion algebras with certain optimal embeddings
# - https://ia.cr/2023/106 Lemma 10

from sage.algebras.quatalg.quaternion_algebra import basis_for_quaternion_lattice
bfql = lambda els: basis_for_quaternion_lattice(els, reverse=True)

Quat1, (i,j,k) = QuaternionAlgebra(-1, -p).objgens()
assert Quat1.discriminant() == p         # ramifies correctly

O0mat = matrix([list(g) for g in [Quat1(1), i, (i+j)/2, (1+k)/2]])
O0 = Quat1.quaternion_order(list(O0mat))

orders = [ (1, identity_matrix(QQ,4), O0mat, i, O0mat, vector((1,0,0,0))) ]

q = ZZ(1)
while len(orders) < num:
    q = next_prime(q)
    if q % 4 != 1:  # restricting to q ≡ 1 (mod 4)
        continue

    Quatq, (ii,jj,kk) = QuaternionAlgebra(-q, -p).objgens()
    if Quatq.discriminant() != p:       # ramifies incorrectly
        continue

    x, y = QuadraticForm(QQ, 2, [1,0,p]).solve(q)
    gamma = x + j*y
    assert gamma.reduced_norm() == q
    ims1 = [Quat1(1), i*gamma, j, k*gamma]
    assert ims1[1]**2 == -q
    assert ims1[2]**2 == -p
    assert ims1[1]*ims1[2] == ims1[3]
    assert ims1[2]*ims1[1] == -ims1[3]
    # (1,ii,jj,kk)->ims1 is an isomorphism Quatq->Quat1
    iso1q = ~matrix(map(list, ims1))

    r = min(map(ZZ, Mod(-p, 4*q).sqrt(all=True)))

    basq = [
            Quatq(1),
            ii,
            (1 + jj) / 2,
            (r + jj) * ii / 2 / q,
        ]

    Oq = Quatq.quaternion_order(basq)
    assert Oq.discriminant() == p   # is maximal

    mat1 = matrix(map(list, basq)) * ~iso1q
    O1 = Quat1.quaternion_order(list(mat1))
    assert O1.discriminant() == p   # is maximal
    assert j in O1                  # p-extremal

    # look for an odd connecting ideal
    I = O0 * O1
    I *= I.norm().denominator()
    assert I.is_integral()
    for v in IntegralLattice(I.gram_matrix()).enumerate_short_vectors():
        elt = sum(c*g for c,g in zip(v,I.basis()))
        if ZZ(elt.reduced_norm() / I.norm()) % 2:
            break
    I = I * (elt.conjugate() / I.norm())
    assert I.is_integral()
    assert I.norm() % 2
    assert I.left_order() == O0

    O1_ = I.right_order()
    assert O1_.unit_ideal() == elt * O1 * ~elt
    idl1 = matrix(map(list, I.basis()))

    # q
    # isomorphism from (-1,-p) algebra to (-q,-p) algebra
    # basis of maximal order O₁ in (-1,-p) algebra
    # element sqrt(-q) in O₁ in (-1,-p) algebra
    # basis of connecting ideal I from O₀ in (-1,-p) algebra
    # element γ such that I has right order γ O₁ γ^-1
    orders.append((q, iso1q, mat1, ims1[1], idl1, vector(elt)))

