from sage.all import ZZ, GF, EllipticCurve, parallel

def even_torsion_basis_E0(E0, f):
    """
    For the case when A = 0 we can't use the entangled basis algorithm
    so we do something "stupid" to simply get something canonical
    """
    assert E0.a_invariants() == (0, 0, 0, 1, 0)

    Fp2 = E0.base_ring()
    p = Fp2.characteristic()

    def points_order_two_f():
        """
        Compute a point P of order 2^f with x(P) = 1 + i*x_im
        """
        x_im = 0
        while True:
            x_im += 1
            x = Fp2([1, x_im])
            if not E0.is_x_coord(x):
                continue
            # compares a+bi <= c+di iff (a,b) <= (c,d) as tuples, where integers
            # modulo p are compared via their minimal non-negative representatives
            P = min(E0.lift_x(x, all=True), key = lambda pt: list(pt.y()))
            P.set_order(multiple=p+1)
            if P.order() % (1 << f) == 0:
                P *= P.order() // (1 << f)
                P.set_order(1 << f)
                yield P

    pts = points_order_two_f()
    P = next(pts)
    for Q in pts:
        # Q is picked to be in E[2^f] AND we must ensure that
        # <P, Q> form a basis, which is the same as e(P, Q) having
        # full order 1 << f.
        e = P.weil_pairing(Q, 1 << f)
        if e ** (1 << f - 1) == -1:
            break

    # Finally we want to make sure Q is above (0, 0)
    P2 = (1 << f - 1) * P
    Q2 = (1 << f - 1) * Q
    if Q2 == E0(0, 0):
        pass
    elif P2 == E0(0, 0):
        P, Q = Q, P
    else:
        Q += P

    assert P.weil_pairing(Q, 1 << f) ** (1 << f - 1) == -1
    assert (1 << f - 1) * Q == E0(0, 0)

    return P, Q


if __name__ == "__main__":
    # p, f = 5 * 2**248 - 1, 248
    # p, f = 65 * 2**376 - 1, 376
    p, f = 27 * 2**500 - 1, 500
    print(f"p = {ZZ(p+1).factor()} - 1")
    Fp2 = GF(p**2, modulus=[1, 0, 1], names="i")
    E = EllipticCurve(Fp2, [1, 0])
    E.set_order((p + 1) ** 2)

    P, Q = even_torsion_basis_E0(E, f)
    print(f"{P = }")
    print(f"{Q = }")

    assert P.order() == 1 << f
    assert Q.order() == 1 << f
    e = P.weil_pairing(Q, 1 << f)
    assert e ** (1 << f - 1) == -1
    print("all good")
