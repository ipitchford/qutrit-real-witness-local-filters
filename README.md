# Two-qutrit entanglement invisible to real witnesses under every local filter

*A counterexample to both parts of the Shen–Chen–Bian conjecture, its smallest-dimension local-unitary case, and a
criterion for every dimension*

Version 0.1.0-candidate. DOI [10.5281/zenodo.24540619](https://doi.org/10.5281/zenodo.24540619) (all versions:
[10.5281/zenodo.24540618](https://doi.org/10.5281/zenodo.24540618)). Release page:
https://evidencepress.org/releases/qutrit-real-witness-local-filters/

**Status: unrefereed research candidate.** Author: Anonymous. Text CC0-1.0; code MIT (see [LICENSE.md](LICENSE.md)).

## What this is about

**The problem.** To show that two quantum particles are entangled, experimenters measure an *entanglement witness*: a
quantity whose average can only go negative for entangled states. A witness is written as a table of numbers (a matrix
in a chosen basis of each particle), and in general those numbers are complex. Real-number tables are simpler, and a
natural hope is that nothing is lost by using only real ones, as long as each laboratory may first pick a good local
basis or apply a *local filter* (an invertible local operation, carried out with some probability of success, that
reshapes its particle's state).

**What was known.** Shen, Chen and Bian (Physical Review A, 2025, Conjecture 9) conjectured exactly this, in two
strengths: (i) a local change of basis suffices (local unitaries); (ii) local filters suffice (SLOCC). Part (i) implies
part (ii). They proved part (i) for every state whose partial transpose has a negative eigenvalue, and both parts for
some two-qutrit families. For a qubit pair and for a qubit with a qutrit (total dimension below eight) the conjecture
holds automatically, because there every entangled state has a partial transpose with a negative eigenvalue.

**What is new.**
1. **Both parts are false** (Theorem 8). An explicit full-rank two-qutrit state (two three-level systems) is detected by
   a complex witness, but every real witness, after any invertible local filters and in any limit of such operations,
   has average above 8.5×10⁻¹⁰ on it. The same holds for every state within Frobenius distance 8×10⁻¹⁰.
2. **Part (i) fails in the smallest possible size** (Theorem 10). For a qubit paired with a four-level system (2×4,
   total dimension eight) there are states that a filtered real witness detects but no locally rotated
   real witness does. Here filters help but basis changes do not. Whether part (ii) also fails for 2×4 is open.
3. **A criterion in every dimension** (Theorem 6). The conjecture holds exactly when every "corner" (exposed point) of
   the set of witnesses can be approximated by filtered real ones. When it fails, it fails for an open set of states.
4. **A test for a known family** (Theorem 9). A degree-nine polynomial decides which members of an exposed family of
   witnesses due to Ha and Kye have filtered real versions.
5. **The fix** (Section 8, Proposition 11). One extra local qubit restores detection; a second copy of the state can
   serve as that qubit when filters are allowed.

**Idea of the proof.** Take a witness that touches the separable states at only one point (an *exposed* witness). If
that witness is not a limit of filtered real ones, nudging the touching separable state slightly towards entanglement
gives states the witness detects but every filtered real witness misses. To show that the chosen witness is not
filtered real, the paper uses a polynomial built from the matrix that changes only by a positive factor under any
local filter and is real for real matrices; for this witness its imaginary part is not zero.

**Limits.** The explicit state is extremely close to separable: the certified ball has radius 8×10⁻¹⁰ and the noise
threshold is about 1.6×10⁻⁷. This is a mathematical statement about which observables suffice, not an experimental
proposal. The witness and the separable state come from earlier work of Ha and Kye. Experiments can still detect such
states with several real measurements and classical post-processing. Further open problems are in Section 9 of the
paper.

## Contents

| File | What it is |
|---|---|
| [manuscript.pdf](manuscript.pdf), [manuscript.tex](manuscript.tex) | The paper (12 pages; the figure is inline TikZ) |
| [state.json](state.json) | The exact state ρ_ε (Gaussian-integer entries over the denominator 57599991) and the witness W |
| [code/verify_witness.py](code/verify_witness.py) | SymPy, exact: rebuilds W and the separable state from the Ha–Kye data, checks the positivity conditions, spectra, the Gram certificate, the invariants and every constant of Section 6, the noise threshold and state.json (120 checks) |
| [code/verify_chirality.py](code/verify_chirality.py) | Standard library, exact: expands the degree-nine invariant and its antisymmetric factorisation |
| [code/independent_tensor_check.py](code/independent_tensor_check.py) | NumPy, exact modular arithmetic: an independent check of the same polynomial identity |
| [code/verify_qubit_ququart.py](code/verify_qubit_ququart.py) | SymPy, exact: rebuilds the 2×4 map from the Ha–Kye formula, checks tr([X,Y]³) = −840i and a real control |
| [code/check_block_positivity.py](code/check_block_positivity.py) | NumPy, floating point: a sanity cross-check of the cited block-positivity theorems, not part of any proof |
| [results/](results/) | Recorded outputs of the five programs, and `verification_witness.json` |
| [requirements.txt](requirements.txt) | sympy==1.14.0, numpy==2.3.5 |
| [refs.bib](refs.bib) | Bibliography (Crossref and arXiv records) |
| [AI_INDEX.md](AI_INDEX.md) | Machine-oriented map of claims, evidence and limits |
| [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md), [RESEARCH_GATES.json](RESEARCH_GATES.json) | Internal review and research-gate records |
| [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), [ENVIRONMENT.txt](ENVIRONMENT.txt) | Release replay and environment |
| [MANIFEST.sha256](MANIFEST.sha256) | SHA-256 of every file |

## Replay

From the package root, with Python 3.12, SymPy 1.14.0 and NumPy 2.3.5:

```sh
OPENBLAS_NUM_THREADS=1 python3 -I -B -O code/verify_witness.py --output results/verification_witness.json
OPENBLAS_NUM_THREADS=1 python3 -I -B -O code/verify_chirality.py
OPENBLAS_NUM_THREADS=1 python3 -I -B -O code/verify_qubit_ququart.py
OPENBLAS_NUM_THREADS=1 python3 -I -B -O code/independent_tensor_check.py
OPENBLAS_NUM_THREADS=1 python3 -I -B -O code/check_block_positivity.py
```

For example `uv run --python 3.12 --with sympy==1.14.0 --with numpy==2.3.5 python3 -I -B -O code/verify_witness.py`.
Each program takes a few seconds and exits with an error on any failed check, also under `-O`. Expected final lines:
`ALL 120 CHECKS PASSED`, `antisymmetric factorization: PASS`, `Commutator trace: -840*I` then `ALL CHECKS PASSED`,
`Integer polynomial identity: PASS`, `ALL CHECKS PASSED (numerical sanity check only)`.

**Comparison policy (declared before tagging; the same in ENVIRONMENT.txt, REPLAY_RECEIPT.md and the CI workflow).**
Every certifying check is exact. The outputs of `verify_chirality.py`, `verify_qubit_ququart.py` and
`independent_tensor_check.py` must be identical line by line. In `verify_witness.py` every line is identical except the
two floating-point display lines (`alpha:` and `lower bound on filtered real expectations:`), whose values must agree
within 1e-9 relative; in `results/verification_witness.json` every field is identical except the float display fields
(`tau_r`, `cubic_change_bound`, `alpha`, `filtered_real_expectation_lower_bound`, `psd_margin`, `ppt_margin`), which must
agree within 1e-9 relative. `check_block_positivity.py` is a floating-point sanity check: its two minima must each be at
least −1e-9 (they are about −1e-15 and depend on the linear-algebra library), and its final line must match exactly.
Interpreter version fields are not compared. CI runs all five programs on Linux under this policy, plus a negative
control.

## Assurance

Written proofs with exact computer checks. Internal AI-assisted review only: a combined five-role Claude review panel,
the lead's proof audit and clean replay, one revision round with extensions, a verification re-review and a proof audit
of the new and rewritten proofs, and a clean release replay with a negative control. No Astra/Codex review took place.
This is not peer review, independent reproduction or formal verification. See [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).

## Cite

Anonymous (2026). *Two-qutrit entanglement invisible to real witnesses under every local filter* (0.1.0-candidate).
Evidence Press. https://doi.org/10.5281/zenodo.24540619
