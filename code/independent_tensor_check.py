#!/usr/bin/env python3
"""Independent dense tensor-network check of the two-qutrit witness invariant.

Requires NumPy. It uses no sparse monomial or coefficient enumeration.
Run with OPENBLAS_NUM_THREADS=1 python3 code/independent_tensor_check.py.
Every check raises SystemExit on failure, so it also fails under python -O.
"""
from itertools import combinations, permutations
from math import prod
import numpy as np

# One-based copy indices. Each consecutive group of three Index symbols
# is contracted against an epsilon tensor in the listed order.
orders = [
    (1,2,3,4,5,6,7,8,9), (1,6,4,9,5,7,3,2,8),
    (2,3,7,6,1,8,9,5,4), (1,4,9,7,5,2,8,3,6),
]
indices = [[(copy, leg) for leg in range(4)] for copy in range(1,10)]
for leg, order in enumerate(orders):
    for start in range(0,9,3):
        indices.append([(copy, leg) for copy in order[start:start+3]])

# Build a pairwise contraction tree that limits intermediate tensor rank.
plan = []
while len(indices) > 1:
    options = []
    for i,j in combinations(range(len(indices)),2):
        shared = set(indices[i]) & set(indices[j])
        if shared:
            result_rank = len(indices[i]) + len(indices[j]) - 2*len(shared)
            options.append((result_rank,i,j))
    _,i,j = min(options)
    right = indices.pop(j)
    left = indices.pop(i)
    shared = sorted(set(left) & set(right))
    plan.append((i,j,[left.index(k) for k in shared],[right.index(k) for k in shared]))
    indices.append([k for k in left+right if k not in shared])

epsilon = np.zeros((3,3,3),dtype=np.int64)
for perm in permutations(range(3)):
    epsilon[perm] = (-1)**sum(x > y for x,y in combinations(perm,2))



def check(condition, label):
    """Explicit check that also runs under python -O."""
    if not bool(condition):
        raise SystemExit(f"FAILED: {label}")

def contraction(a,b,c,u,v,prime):
    w = np.zeros((3,3,3,3),dtype=np.int64)
    for i,j in np.ndindex(3,3):
        w[i,j,i,j] = (a,b,c)[(j-i)%3]
        if i != j:
            w[i,i,j,j] = u if (j-i)%3==1 else v
    tensors = [w%prime]*9 + [epsilon%prime]*12
    for i,j,left_axes,right_axes in plan:
        right = tensors.pop(j)
        left = tensors.pop(i)
        tensors.append(np.tensordot(left,right,axes=(left_axes,right_axes))%prime)
    return int(tensors[0])


def reconstruct(primes,values):
    modulus = prod(primes)
    answer = sum(x*(modulus//p)*pow(modulus//p,-1,p) for p,x in zip(primes,values))%modulus
    return answer if 2*answer < modulus else answer-modulus


def main():
    # Separate real and imaginary components using both embeddings of Z[i].
    primes = [101,241,1009,1013,1021,1033,1049]
    real,imag = [],[]
    for prime in primes:
        root = next(a for a in range(prime) if (a*a+1)%prime==0)
        plus,minus = [contraction(3,1,4,-4-3*h,-4+3*h,prime) for h in (root,-root)]
        real.append((plus+minus)*pow(2,-1,prime)%prime)
        imag.append((plus-minus)*pow(2*root,-1,prime)%prime)
    # Each epsilon has six possible nonzero index configurations and |W_ijkl|≤5.
    check(prod(primes) > 2*(6**12)*(5**9), "independent_tensor_check.py check 1")
    value = reconstruct(primes,real),reconstruct(primes,imag)
    check(value==(188339589,3582306), "independent_tensor_check.py check 2")
    print('Exact invariant:',value)

    # For v=1 both sides have degree at most 9 in the other four variables.
    # The 715 nonnegative integer lattice points with sum at most 9 determine
    # this polynomial over all four fields (Newton interpolation). Homogeneity
    # recovers the identity in v. Total coefficient difference is bounded by
    # two copies of 6^12 plus the l1 coefficient norm bound 264 of the product.
    for prime in primes[:4]:
        for i,j,k,l in combinations(range(13),4):
            a,b,c,u,v=i,j-i-1,k-j-1,l-k-1,1
            factor = sum(x**3 for x in (a,b,c,u,v)) - 3*a*b*c - 3*a*u*v
            target = -6*(b**3-c**3)*(u**3-v**3)*factor
            diff = contraction(a,b,c,u,v,prime)-contraction(a,b,c,v,u,prime)
            check((diff-target)%prime==0, "independent_tensor_check.py check 3")
        print('Polynomial identity on all 715 interpolation points modulo',prime)
    check(prod(primes[:4]) > 2*(6**12)+264, "independent_tensor_check.py check 4")
    print('Integer polynomial identity: PASS')


if __name__=='__main__':
    main()
