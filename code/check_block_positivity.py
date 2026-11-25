"""Floating-point sanity check (not a proof) of block positivity of W and of the 2x4 matrix F.

Block positivity of W is a theorem of Ha and Kye (2012, Theorem 2.2); that of F follows from
Ha and Kye (2016, Theorem 3.2). This program only cross-checks those theorems numerically:
for x on the unit sphere of the first factor, the compressed matrix
    W_x = (x^dagger (x) I) W (x (x) I)
must be positive semidefinite. It samples 200000 random unit x (fixed seed) plus the 16+12
first-factor vectors of the exposing family and reports the smallest eigenvalue found.
Declared tolerance: the minimum must be >= -1e-9 (eigenvalues are O(1)), on any platform.
"""
import numpy as np


def check(condition, label):
    if not bool(condition):
        raise SystemExit(f"FAILED: {label}")


def witness():
    w = np.zeros((9, 9), dtype=complex)
    for i in range(3):
        for j in range(3):
            w[3*i+j, 3*i+j] = (3, 1, 4)[(j - i) % 3]
    for i in range(3):
        j = (i + 1) % 3
        w[4*i, 4*j] = -4 - 3j
        w[4*j, 4*i] = -4 + 3j
    return w


def woronowicz_choi():
    def phi(x):
        xx, y, z, w = x[0, 0], x[0, 1], x[1, 0], x[1, 1]
        return np.array([
            [4*xx - 2*(y + z) + 3*w, -2*xx + 2*z, 0, 0],
            [-2*xx + 2*y, 2*xx, z, 0],
            [0, y, 2*w, -2*z - w],
            [0, 0, -2*y - w, 4*xx + 2*w]], dtype=complex)
    out = np.zeros((8, 8), dtype=complex)
    for i in range(2):
        for j in range(2):
            e = np.zeros((2, 2))
            e[i, j] = 1
            out[4*i:4*i+4, 4*j:4*j+4] = phi(e)
    return out


def min_compressed_eig(v, m, n, samples, rng, extra=()):
    blocks = v.reshape(m, n, m, n)
    xs = rng.normal(size=(samples, m)) + 1j*rng.normal(size=(samples, m))
    xs = np.vstack([xs] + [np.asarray(e, dtype=complex)[None, :] for e in extra])
    xs /= np.linalg.norm(xs, axis=1, keepdims=True)
    comp = np.einsum('si,ijkl,sk->sjl', xs.conj(), blocks, xs)
    return float(np.linalg.eigvalsh(comp).min())


def main():
    rng = np.random.default_rng(20261124)
    phases = (1, 1j, -1, -1j)
    extra = [(1, u, v) for u in phases for v in phases]
    for i in range(3):
        for u in phases:
            x = [0, 0, 0]
            x[i], x[(i + 1) % 3] = (1 + 1j)*u, 1
            extra.append(tuple(x))
    w_min = min_compressed_eig(witness(), 3, 3, 200000, rng, extra)
    f_min = min_compressed_eig(woronowicz_choi(), 2, 4, 200000, rng)
    print(f"min eigenvalue of compressed W: {w_min:.3e}")
    print(f"min eigenvalue of compressed F: {f_min:.3e}")
    check(w_min >= -1e-9, "W block-positive on samples")
    check(f_min >= -1e-9, "F block-positive on samples")
    check(abs(w_min) < 1e-9, "W has product zeros (exposing family attains 0)")
    print("ALL CHECKS PASSED (numerical sanity check only)")


if __name__ == "__main__":
    main()
