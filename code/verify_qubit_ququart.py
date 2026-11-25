"""Local-unitary obstruction in the least possible total dimension (two by four).

Ha and Kye, arXiv:1410.5545, Theorem 3.2, prove that the map Phi[a,b,c,d] from M2 to M4 is
exposed and indecomposable whenever a,b,c,d > 0 and ab > 1. At (a,b,c,d) = (2,2,2,1) it is,
as they note, Woronowicz's example (with 3d in the (1,1) entry). This program rebuilds the map
from their general formula, forms the Choi matrix F, applies the complex local shear
B = I + i E_01 on the four-dimensional factor, and evaluates the local-unitary invariant
tr([X,Y]^3) with X = tr_A F', Y = tr_A F'^2.
Every check raises SystemExit on failure, so the program also fails under python -O.
"""

from sympy import I, Matrix, Rational, expand, eye, kronecker_product, sqrt, trace, zeros


def check(condition, label):
    """Explicit check that also runs under python -O."""
    if not bool(condition):
        raise SystemExit(f"FAILED: {label}")


def ha_kye_map(a, b, c, d):
    """Phi[a,b,c,d] of Ha-Kye (2016), with e,f,g,h,k fixed by their relations (1)."""
    # (ab-1)(e, f) = a(c+d)(c, d), g^2 = acd, h = be - c^2, k = bf - d^2.
    e = Rational(a*(c + d)*c, a*b - 1)
    f = Rational(a*(c + d)*d, a*b - 1)
    g = sqrt(a*c*d)
    h = b*e - c**2
    k = b*f - d**2

    def phi(x):
        xx, y, z, w = x[0, 0], x[0, 1], x[1, 0], x[1, 1]
        return Matrix([
            [h*xx - c*d*(y + z) + k*w, -g*xx + g*z, 0, 0],
            [-g*xx + g*y, a*xx, z, 0],
            [0, y, b*w, -c*z - d*w],
            [0, 0, -c*y - d*w, e*xx + f*w],
        ])
    return phi


def phi_displayed(x):
    """The map as displayed in the manuscript, Section 7."""
    xx, y, z, w = x[0, 0], x[0, 1], x[1, 0], x[1, 1]
    return Matrix([
        [4*xx - 2*(y + z) + 3*w, -2*xx + 2*z, 0, 0],
        [-2*xx + 2*y, 2*xx, z, 0],
        [0, y, 2*w, -2*z - w],
        [0, 0, -2*y - w, 4*xx + 2*w],
    ])


def unit(i, j):
    x = zeros(2)
    x[i, j] = 1
    return x


def choi(phi):
    return Matrix.vstack(*[Matrix.hstack(*[phi(unit(i, j)) for j in range(2)]) for i in range(2)])


def marginal(v):
    return (v[:4, :4] + v[4:, 4:]).applyfunc(expand)


def main():
    a, b, c, d = 2, 2, 2, 1
    check(a*b > 1, "Ha-Kye (2016) Theorem 3.2 hypothesis ab > 1")
    w = choi(ha_kye_map(a, b, c, d))
    check(w == choi(phi_displayed), "displayed map equals Ha-Kye Phi[2,2,2,1]")
    check(w == w.H and w == w.conjugate(), "Choi matrix F is real symmetric")
    shear = eye(4)
    shear[0, 1] = I
    k = kronecker_product(eye(2), shear)
    v = (k*w*k.H).applyfunc(expand)
    x, y = marginal(v), marginal(v**2)
    comm = x*y - y*x
    value = expand(trace(comm**3))
    check(value == -840*I, "tr([X,Y]^3) = -840 i for F' = (I x B) F (I x B)^dagger")
    # Control: the unsheared real F gives zero.
    x0, y0 = marginal(w), marginal(w**2)
    check(expand(trace((x0*y0 - y0*x0)**3)) == 0, "control: invariant vanishes on the real F")
    print("Commutator trace:", value)
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
