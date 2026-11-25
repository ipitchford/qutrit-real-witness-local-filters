"""Exact certificate for a two-qutrit state missed by every locally filtered real witness.

Run from the package root with Python 3 and SymPy:
    python3 code/verify_witness.py [--output results/verification_witness.json]

Indices i,j label |i,j> in lexicographic order (first index changes slowly).
Hermitian conjugation is .H. States and witnesses are paired by the ordinary trace.
Every check raises SystemExit on failure, so the program also fails under python -O.
"""

import json
import sys
from itertools import product
from pathlib import Path

from sympy import (I, Matrix, Rational, SparseMatrix, conjugate, expand, eye, im,
                   kronecker_product, nsimplify, re, simplify, sqrt, trace, zeros)

CHECKS = []


def check(condition, label):
    if not bool(condition):
        raise SystemExit(f"FAILED: {label}")
    CHECKS.append(label)


def flat(i, j):
    return 3 * i + j


def T(a, b, c, u):
    """Diagonal (a,b,c)[(j-i) mod 3] at |ij>; u from row |ii> to column |i+1,i+1>."""
    m = zeros(9)
    for i, j in product(range(3), repeat=2):
        m[flat(i, j), flat(i, j)] = (a, b, c)[(j - i) % 3]
    for i in range(3):
        j = (i + 1) % 3
        m[flat(i, i), flat(j, j)] = u
        m[flat(j, j), flat(i, i)] = conjugate(u)
    return m


def ha_kye_choi(a, b, c, phase):
    """Choi matrix sum_ij e_ij (x) Phi[a,b,c;theta](e_ij) of the map in Ha-Kye (2013), Sec. 1.

    phase = e^{i theta}. Phi(X) has diagonal (a x11+b x22+c x33, c x11+a x22+b x33,
    b x11+c x22+a x33) and off-diagonal -e^{i theta} x12, -e^{-i theta} x13,
    -e^{-i theta} x21, -e^{i theta} x23, -e^{i theta} x31, -e^{-i theta} x32.
    """
    ph, cph = phase, conjugate(phase)

    def phi(x):
        return Matrix([
            [a*x[0, 0] + b*x[1, 1] + c*x[2, 2], -ph*x[0, 1], -cph*x[0, 2]],
            [-cph*x[1, 0], c*x[0, 0] + a*x[1, 1] + b*x[2, 2], -ph*x[1, 2]],
            [-ph*x[2, 0], -cph*x[2, 1], b*x[0, 0] + c*x[1, 1] + a*x[2, 2]],
        ])
    out = zeros(9)
    for i, j in product(range(3), repeat=2):
        e = zeros(3)
        e[i, j] = 1
        out[3*i:3*i+3, 3*j:3*j+3] = phi(e)
    return out


def ha_kye_rho0(t, phase):
    """Unnormalised separable state rho_0(theta,t) of Ha-Kye (2013), Section 4, built from
    their product vectors z1(1,i^k,i^l) (weight 1/16) and z_k(e^{2 pi i l/3}) (weight 1/(3t^2))."""
    omega = Rational(-1, 2) + sqrt(3)*I/2
    st = sqrt(t)
    rho = zeros(9)
    for k, l in product(range(1, 5), repeat=2):
        x = Matrix([1, I**k, I**l])
        z = kronecker_product(x, x.conjugate())
        rho += z * z.H / 16
    for l in range(1, 4):
        beta = omega**l
        bb = conjugate(beta)
        vecs = [
            (Matrix([0, st*beta, 1]), Matrix([0, st*bb, t*phase])),
            (Matrix([1, 0, st*beta]), Matrix([t*phase, 0, st*bb])),
            (Matrix([st*beta, 1, 0]), Matrix([st*bb, t*phase, 0])),
        ]
        for x, y in vecs:
            z = kronecker_product(x, y)
            rho += z * z.H / (3*t**2)
    return rho.applyfunc(lambda v: nsimplify(expand(v)))


def zero_vectors():
    phases = (1, I, -1, -I)
    for u, v in product(phases, repeat=2):
        x = Matrix([1, u, v])
        yield x, x.conjugate()
    for i, u in product(range(3), phases):
        j = (i + 1) % 3
        x, y = zeros(3, 1), zeros(3, 1)
        x[i], x[j] = (1 + I) * u, 1
        y[i], y[j] = (1 - I) * conjugate(u), 2 * (4 - 3*I)/5
        yield x, y


def exposing_state_and_constraints():
    d = zeros(9)
    rows = []
    norms = []
    for x, y in zero_vectors():
        z = kronecker_product(x, y)
        d += z * z.H
        norms.append(expand((x.H*x)[0] + (y.H*y)[0]))
        for k in range(3):
            e = eye(3)[:, k]
            for v in (kronecker_product(e, y), kronecker_product(x, e)):
                rows.append(kronecker_product(v.H, z.T).applyfunc(expand))
    return d.applyfunc(expand), SparseMatrix(Matrix.vstack(*rows)), norms


def partial_transpose(h):
    return Matrix(9, 9, lambda r, s: h[flat(r // 3, s % 3), flat(s // 3, r % 3)])


def covariant(x, dual):
    """C(X,Y)_{jb,ck} = sum_{i,l} X_{ij,kl} Y_{bl,ic}  (manuscript, Section 4)."""
    return Matrix(9, 9, lambda r, s: expand(sum(
        x[flat(i, r//3), flat(s % 3, l)] * dual[flat(r % 3, l), flat(i, s//3)]
        for i, l in product(range(3), repeat=2))))


def verify_gram(q):
    identity_vector = SparseMatrix(81, 1, {(k*9+k, 0): 1 for k in range(9)})
    g = (q.H * q + identity_vector * identity_vector.T - Rational(2, 5)*eye(81)).applyfunc(expand)
    check(g == g.H, "Gram matrix Hermitian")
    components = g.connected_components()
    minima = []
    for indices in components:
        a = g.extract(indices, indices)
        while a.rows > 0:
            pivot = a[0, 0]
            check(pivot.is_Rational and pivot > 0, "positive rational Schur pivot")
            minima.append(pivot)
            a = (a[1:, 1:] - a[1:, 0] * a[0, 1:] / pivot).applyfunc(expand)
    return sorted(map(len, components)), min(minima)


def main():
    out_path = None
    if "--output" in sys.argv:
        out_path = Path(sys.argv[sys.argv.index("--output") + 1])
    report = {}

    # 1. The witness W and the Ha-Kye data at e^{i theta} = (4+3i)/5, t = 2.
    phase = (4 + 3*I)/5
    t = 2
    w = T(3, 1, 4, -4 - 3*I)
    check(w == w.H, "W Hermitian")
    check(re(phase) > Rational(1, 2), "cos(theta) = 4/5 > 1/2, so |theta| < pi/3 and p_theta = 2 cos(theta)")
    p_theta = 2*re(phase)
    s = 1 - t + t**2
    a_hk, b_hk, c_hk = 1 - (p_theta-1)*t/s, (p_theta-1)*t**2/s, (p_theta-1)/s
    check((a_hk, b_hk, c_hk) == (Rational(3, 5), Rational(4, 5), Rational(1, 5)), "Ha-Kye a(2), b(2), c(2)")
    check(a_hk + b_hk + c_hk == p_theta and b_hk*c_hk == (1 - a_hk)**2,
          "Ha-Kye (2012) Theorem 2.2 positivity conditions hold (with equality)")
    check(w == 5*ha_kye_choi(a_hk, b_hk, c_hk, phase), "W = 5 x Choi matrix of Phi_theta(2)")
    check((w.det(), trace(w), trace(w**2)) == (-7040, 24, 228), "det, trace and trace of square of W")
    check(w.eigenvals() == {-5: 1, 1: 3, 4: 3, 7 - 3*sqrt(3): 1, 7 + 3*sqrt(3): 1}, "spectrum of W")

    # 2. The exposing separable operator D.
    d, q, norms = exposing_state_and_constraints()
    check(d == T(48, 48, 24, Rational(144, 5) + Rational(48, 5)*I), "28-term sum D equals T(48,48,24;(144+48i)/5)")
    check(d == (16*ha_kye_rho0(t, phase).T).applyfunc(expand), "D = 16 x transpose of Ha-Kye rho_0(theta,2)")
    check(trace(d) == 360 and expand(trace(w*d)) == 0, "tr D = 360 and tr(WD) = 0")
    check(max(norms) == 9, "max of ||x||^2+||y||^2 over the 28 pairs is 9")
    check(q.shape == (168, 81) and (q*Matrix(list(w))).applyfunc(expand) == zeros(168, 1),
          "168 first-order rows annihilate vec(W)")
    lam_d = min(d.eigenvals())
    lam_dg = min(partial_transpose(d).eigenvals())
    check(lam_d == Rational(96, 5) - 48*sqrt(3)/5 and lam_d > 2, "lambda_min(D) = (96-48 sqrt3)/5 > 2")
    check(lam_dg == 36 - 12*sqrt(185)/5 and lam_dg > 3, "lambda_min(D^Gamma) = 36-12 sqrt185/5 > 3")
    blocks, smallest_pivot = verify_gram(q)
    check(blocks == [1]*6 + [2]*6 + [8]*6 + [15], "Gram block sizes 1,2,8,15 with multiplicities 6,6,6,1")
    check(smallest_pivot == Rational(6849698078814, 29474774167205), "smallest Schur pivot")
    report["min_eig_D"] = str(lam_d)
    report["min_eig_D_Gamma"] = str(lam_dg)
    report["gram_blocks"] = blocks
    report["smallest_schur_pivot"] = str(smallest_pivot)

    # 3. The filter invariant chi.
    m = covariant(w, w.adjugate())
    invariant = expand(trace(m**3))
    check(invariant == 5007168700416 + 151720116387840*I, "tr M(W)^3")
    h = w / 24
    hinv = h.inv()
    m0 = covariant(h, hinv)
    im_l3 = simplify(im(expand(trace(m0**3))))
    check(im_l3 == Rational(-84184191, 193600), "Im tr L(H)^3 = -84184191/193600 (computed)")
    check(simplify(im(invariant)/w.det()**3) == im_l3, "Im tr L(H)^3 = Im tr M(W)^3 / det(W)^3")
    check(abs(im_l3) > 434, "|Im tr L(H)^3| > 434")
    zc = 33024 + 17280*I
    check(m[0, 4] == zc and m[4, 8] == zc and m[8, 0] == zc and m[4, 0] == conjugate(zc)/4
          and m[0, 0] == -3072, "circulant block of M(W): diagonal -3072, z forward, conj(z)/4 reverse")
    check(expand(3*(1 - Rational(1, 64))*im(expand(zc**3))) == 151720116387840, "chi(W) = 3(1-4^-3) Im z^3")
    report["trace_M_W_cubed"] = str(invariant)
    report["Im_trace_L_H_cubed"] = str(im_l3)

    # 4. Norm constants of Section 4.
    check(max(abs(v) for v in hinv.eigenvals()) == 24, "||H^{-1}||_op = 24")
    check(expand(trace(m0*m0.H)) == Rational(5157567, 48400) < 11**2, "||L(H)||_F^2 < 11^2")
    check(expand(trace((m0**2)*(m0**2).H)) == Rational(7157861244459, 2342560000) < 56**2, "||L(H)^2||_F^2 < 56^2")
    h_pt_norm = max(abs(v) for v in partial_transpose(h).eigenvals())
    hinv_pt_norm = max(abs(v) for v in partial_transpose(hinv).eigenvals())
    check(h_pt_norm == (5 + sqrt(109))/48 and h_pt_norm < Rational(322, 1000), "||H^Gamma||_op < 0.322")
    check(hinv_pt_norm == 15 + 3*sqrt(36601)/55 and hinv_pt_norm < 26, "||(H^{-1})^Gamma||_op < 26")
    radius = Rational(1, 140)
    tau = (26 + (Rational(322, 1000)+radius)*24**2/(1-24*radius))*radius
    check(tau < Rational(1824, 1000), "tau_r < 1.824")
    change = 3*56*tau + 3*11*tau**2 + tau**3
    check(change < 424, "cubic-trace change < 424")
    eps = Rational(1, 160000)
    alpha = Rational(2, 5)*radius**2/(9*(Rational(1, 3)+sqrt(3)*radius))
    check(alpha > Rational(65592, 10**10) and alpha > eps, "alpha > 6.5592e-6 > epsilon")
    lower = (alpha - eps)/(360 - 9*eps)
    check(lower > Rational(85, 10**11), "(alpha-eps)/(360-9 eps) > 8.5e-10")
    report["tau_r"] = float(tau)
    report["cubic_change_bound"] = float(change)
    report["alpha"] = float(alpha)
    report["filtered_real_expectation_lower_bound"] = float(lower)

    # 5. The state rho_eps and its ball of Frobenius radius 8e-10.
    rho = (d - eps*eye(9))/(360 - 9*eps)
    w_exp = expand(trace(w*rho))
    check(trace(rho) == 1 and w_exp == Rational(-8, 19199997), "tr(W rho_eps) = -8/19199997")
    ball = Rational(8, 10**10)
    margin_psd = (lam_d - eps)/(360 - 9*eps)
    margin_ppt = (lam_dg - eps)/(360 - 9*eps)
    check(margin_psd > ball and margin_ppt > ball, "eigenvalue margins of rho_eps and rho_eps^Gamma exceed 8e-10")
    check(eps/(360 - 9*eps) > ball and lower > ball, "a Frobenius perturbation below 8e-10 preserves both witness signs")
    report["tr_W_rho_eps"] = str(w_exp)
    report["psd_margin"] = float(margin_psd)
    report["ppt_margin"] = float(margin_ppt)

    # 6. Noise threshold along rho(p) = (1-p) rho_eps + p I/9.
    p_star = eps/40
    check(((1 - p_star)*rho + p_star*eye(9)/9).applyfunc(expand) == (d/360).applyfunc(expand), "rho(eps/40) = D/360")
    check(expand((1 - p_star)*trace(w*rho) + p_star*trace(w)/9) == 0, "tr(W rho(p*)) = 0")

    # 7. state.json agrees with rho_eps and W.
    root = Path(__file__).resolve().parent.parent
    data = json.loads((root / "state.json").read_text())
    den = data["state_denominator"]
    num = Matrix(9, 9, lambda r, c: data["state_numerator"][r][c][0] + I*data["state_numerator"][r][c][1])
    check(den == 57599991 and (num/den - rho).applyfunc(expand) == zeros(9), "state.json density matrix = rho_eps")
    wj = Matrix(9, 9, lambda r, c: data["witness"][r][c][0] + I*data["witness"][r][c][1])
    check(wj == w, "state.json witness = W")

    print("min eig(D), min eig(D^Gamma):", lam_d, lam_dg)
    print("Gram block sizes:", blocks)
    print("Smallest Schur pivot:", smallest_pivot)
    print("trace M(W,adj W)^3:", invariant)
    print("Im tr L(H)^3:", im_l3)
    print("alpha:", float(alpha))
    print("lower bound on filtered real expectations:", float(lower))
    print("Witness expectation of the state:", w_exp)
    print(f"ALL {len(CHECKS)} CHECKS PASSED")
    report["checks"] = CHECKS
    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=1) + "\n")


if __name__ == "__main__":
    main()
