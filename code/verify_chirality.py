#!/usr/bin/env python3
"""Exact relative-invariant certificate for the two-qutrit witness.

Uses only the Python standard library.  The variables in the polynomial are
the diagonal weights a,b,c (column minus row modulo 3 on a product ket)
and the coherence weights u,v in the forward/reverse orientation.
"""
from collections import Counter
from itertools import product

# Each row gives three ordered triples of copy labels for the epsilon tensors.
TRIPLES = (
    ((0, 1, 2), (3, 4, 5), (6, 7, 8)),  # first ket index
    ((0, 5, 3), (8, 4, 6), (2, 1, 7)),  # second ket index
    ((1, 2, 6), (5, 0, 7), (8, 4, 3)),  # first bra index
    ((0, 3, 8), (6, 4, 1), (7, 2, 5)),  # second bra index
)
# Tuples (ket1,ket2,bra1,bra2), with the index of the associated weight.
SUPPORT = [((i, j, i, j), (j-i) % 3) for i, j in product(range(3), repeat=2)]
SUPPORT += [((i, i, j, j), 2+(j-i) % 3)
            for i, j in product(range(3), repeat=2) if i != j]



def check(condition, label):
    """Explicit check that also runs under python -O."""
    if not bool(condition):
        raise SystemExit(f"FAILED: {label}")

def contraction():
    # Assign all nine factors one by one. Repeated labels within any
    # epsilon triple give zero, so those assignments are discarded at once.
    assignment = [[-1] * 9 for _ in range(4)]
    peers = [[next(t for t in row if slot in t) for slot in range(9)]
             for row in TRIPLES]
    powers = [0] * 5
    coeffs = Counter()
    leaves = [0, 0]  # unsigned counts of negative and positive epsilon products

    def visit(slot):
        if slot == 9:
            inversions = sum(assignment[k][x] > assignment[k][y]
                             for k, row in enumerate(TRIPLES)
                             for t in row for x, y in ((t[0], t[1]), (t[0], t[2]), (t[1], t[2])))
            sign = (-1) ** inversions
            leaves[int(sign > 0)] += 1
            coeffs[tuple(powers)] += sign
            return
        for entry, symbol in SUPPORT:
            if any(entry[k] in (assignment[k][x] for x in peers[k][slot])
                   for k in range(4)):
                continue
            for k in range(4):
                assignment[k][slot] = entry[k]
            powers[symbol] += 1
            visit(slot+1)
            powers[symbol] -= 1
            for k in range(4):
                assignment[k][slot] = -1

    visit(0)
    return Counter({key: c for key, c in coeffs.items() if c}), leaves


def multiply(left, right):
    result = Counter()
    for k, a in left.items():
        for l, b in right.items():
            result[tuple(i+j for i, j in zip(k, l))] += a*b
    return Counter({key: c for key, c in result.items() if c})


def evaluate(coeffs):
    # Gaussian integers, represented as (real, imaginary).
    values = ((3, 0), (1, 0), (4, 0), (-4, -3), (-4, 3))
    real = imag = 0
    for exponent, coefficient in coeffs.items():
        x, y = coefficient, 0
        for n, (r, s) in zip(exponent, values):
            for _ in range(n):
                x, y = x*r-y*s, x*s+y*r
        real += x
        imag += y
    return real, imag


def main():
    coeffs, leaves = contraction()
    delta = Counter(coeffs)
    for (a, b, c, u, v), n in coeffs.items():
        delta[a, b, c, v, u] -= n
    # Verify P(W)-P(conjugate W) as an identity over the polynomial ring:
    # -6(b^3-c^3)(u^3-v^3)(a^3+b^3+c^3+u^3+v^3-3abc-3auv).
    cubed = [tuple(3*int(j == k) for j in range(5)) for k in range(5)]
    first = {cubed[1]: -6, cubed[2]: 6}
    second = {cubed[3]: 1, cubed[4]: -1}
    third = dict.fromkeys(cubed, 1)
    third[(1, 1, 1, 0, 0)] = third[(1, 0, 0, 1, 1)] = -3
    check({k: v for k, v in delta.items() if v} == multiply(multiply(first, second), third), "verify_chirality.py check 1")
    check(evaluate(coeffs) == (188339589, 3582306), "verify_chirality.py check 2")
    print("negative/positive epsilon assignments:", *leaves)
    print("nonzero monomial coefficients:", len(coeffs))
    print("exact P(W): 188339589 + 3582306 i")
    print("antisymmetric factorization: PASS")


if __name__ == "__main__":
    main()
