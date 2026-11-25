# Replay receipt (release 0.1.0-candidate, 25 November 2026)

The publisher replayed the approved candidate from fresh copies before adding the DOI line. Only `manuscript.tex`
(the DOI line under the author), the rebuilt `manuscript.pdf` and the release records changed afterwards; the five
programs, `state.json`, `requirements.txt` and the recorded outputs are byte-identical to the replayed copies (checked
by SHA-256 in the frozen-archive gate).

The interpreter was printed inside the run (see [ENVIRONMENT.txt](ENVIRONMENT.txt)): CPython 3.12.11, SymPy 1.14.0,
NumPy 2.3.5.

## Runs (fresh copies, cwd = package root, macOS arm64, PYTHONDONTWRITEBYTECODE=1, OPENBLAS_NUM_THREADS=1)

Each program ran with `python3 -I -B -O` and again with `python3 -I -B` (no `-O`).

| Program | Exit (both modes) | Time | Output |
|---|---|---|---|
| `code/verify_witness.py --output results/verification_witness.json` | 0 | about 1 s | `ALL 120 CHECKS PASSED`; stdout and JSON byte-identical to the recorded files |
| `code/verify_chirality.py` | 0 | about 2 s | `antisymmetric factorization: PASS`; byte-identical |
| `code/verify_qubit_ququart.py` | 0 | under 1 s | `Commutator trace: -840*I`, `ALL CHECKS PASSED`; byte-identical |
| `code/independent_tensor_check.py` | 0 | about 4 s | `Integer polynomial identity: PASS`; byte-identical |
| `code/check_block_positivity.py` | 0 | about 1 s | minima -1.172e-15 (W) and -2.592e-15 (F), tolerance -1e-9; byte-identical |

Every run wrote nothing to stderr.

## Negative control (fresh copy, `python3 -I -B -O`)

| Mutation | Result |
|---|---|
| `verify_witness.py`: ε changed from 1/160000 to 1/100000 | `FAILED: alpha > 6.5592e-6 > epsilon`, exit 1 |

## Comparison policy (declared before tagging; also in README.md, ENVIRONMENT.txt and the CI workflow)

- `verify_chirality`, `verify_qubit_ququart`, `independent_tensor_check`: identical line by line.
- `verify_witness.txt`: identical except the float display lines `alpha:` and `lower bound on filtered real
  expectations:` (values within 1e-9 relative).
- `verification_witness.json`: identical except the float display fields `tau_r`, `cubic_change_bound`, `alpha`,
  `filtered_real_expectation_lower_bound`, `psd_margin`, `ppt_margin` (within 1e-9 relative).
- `check_block_positivity.txt`: identical outside the numbers; each minimum at least -1e-9; final line identical.
- Interpreter versions are not compared.

The comparison script is `.github/ci/compare_output.py`; it was tested against the release outputs (accepted), against
a copy with a minimum of -2e-9 (rejected) and against a copy with one exact line changed by one unit and a JSON float
changed in the fifth significant digit (rejected).
