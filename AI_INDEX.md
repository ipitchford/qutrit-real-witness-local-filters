# AI index: Two-qutrit entanglement invisible to real witnesses under every local filter (v0.1.0-candidate)

Machine-oriented map of claims, evidence and limits. Human readers should start with [README.md](README.md) ("What
this is about") and the paper [manuscript.pdf](manuscript.pdf) (source [manuscript.tex](manuscript.tex); Section 1 is
a plain introduction with a roadmap). DOI 10.5281/zenodo.24540619 (concept 10.5281/zenodo.24540618).
Status: unrefereed candidate; internal AI-assisted review only.

## Normalisation

Bipartite C^m ⊗ C^n, states and witnesses paired by tr(Vρ). B = trace-one block-positive operators (⟨x⊗y|V|x⊗y⟩ >= 0
for all product vectors). A witness is a block-positive V that is not positive semidefinite. Local filter K = A ⊗ B with
A, B invertible. An operator is filtered real if it equals K S K^† with K a local filter and S real symmetric in the
product basis (Definition 3); R = closure of the filtered real members of B; R_U the same with A, B unitary. Shen-Chen-Bian (SCB, Phys. Rev. A 112, 022415, 2025,
Definition 1): E^LU = {(U⊗V) S (U⊗V)^† : S a real witness, U, V unitary}; E^SLOCC with U, V invertible. SCB
Conjecture 9: (i) every entangled state is detected by some witness in E^LU; (ii) by some witness in E^SLOCC. (i)
implies (ii). The Ha-Kye witness W = 5 × Choi matrix of Phi_theta(2), e^{i theta} = (4+3i)/5; ρ_ε with ε = 1/160000
(state.json, denominator 57599991).

## Claims and where they are proved

| ID | Claim | Location | Status / credit |
|---|---|---|---|
| A (main) | SCB Conjecture 9 fails in both parts on C^3 ⊗ C^3: ρ_ε is full-rank, PPT (positive-definite partial transpose), entangled, tr(W ρ_ε) = -8/19199997, and every element of R has expectation > 8.5×10^-10 on ρ_ε; the same on the Frobenius ball of radius 8×10^-10 | Theorem 8; Sections 4 and 6 | written proof with exact constants; new. W and the separable state D are Ha and Kye's (HK 2012, HK 2013) |
| B | Conjecture 9(i) fails on C^2 ⊗ C^4 (total dimension 8, the minimum): an exposed point of B that is filtered real but outside R_U, so full-rank PPT states are detected in E^SLOCC but not in E^LU | Theorem 10; Section 7 | written proof; tr([X,Y]^3) = -840i checked exactly; new. Whether (ii) fails on 2×4 is open |
| C | Criterion in every dimension: (ii) holds on C^m ⊗ C^n iff every exposed (equivalently extreme) point of B lies in R; same for (i) with R_U; failures are robust (open set of full-rank PPT states) | Theorem 6 (via Lemmas 4, 5); Section 3 | written proof; self-contained farthest-point argument |
| D | Lemma: χ (conjugation-odd relative invariant) vanishes on R | Lemma 7; Section 4 | written proof |
| E | Degree-nine invariant decides which exposed Ha-Kye witnesses T(a, β, γ; u) have filtered real representatives | Theorem 9; Section 5 | written proof; polynomial identity checked exactly (two independent programs) |
| F | One local ancilla qubit restores detection within E^LU; two copies suffice with local filters (second copy acts as a local ancilla); finitely many copies suffice in a fixed basis | Section 8, Proposition 11 | written proof; realification as in McKague-Mosca-Gisin and Renou et al.; LU-only copy number open |
| G | Failures are dimension-specific: an embedded two-qutrit counterexample is detected by an LU-real witness once a local dimension is at least 6 | Section 9 remark | written argument |

## Evidence

| Check | Script | Output | Scope |
|---|---|---|---|
| W from HK data; HK 2012 Thm 2.2 conditions (with equality); D = 16 ρ_0^T and the 28-vector sum; spectra; λ_min(D) = (96-48√3)/5, λ_min(D^Γ) = 36-12√185/5; Gram certificate (block sizes, smallest Schur pivot 6849698078814/29474774167205); tr M(W)^3 = 5007168700416+151720116387840i; Im tr L(H)^3 = -84184191/193600; α > 6.5592e-6 > ε; lower bound 8.5909e-10 > 8.5e-10; ball 8e-10 (PSD margin 7.15e-3, PPT margin 9.32e-3); η* = ε/40; tr(W ρ_ε) = -8/19199997; state.json | [code/verify_witness.py](code/verify_witness.py) | [results/verify_witness.txt](results/verify_witness.txt), [results/verification_witness.json](results/verification_witness.json) | exact (SymPy rationals and algebraic numbers); 120 checks |
| Degree-nine invariant Pi(W) = 188339589 + 3582306 i; 6666/4668 assignments; 84 coefficients; antisymmetric factorisation | [code/verify_chirality.py](code/verify_chirality.py) | [results/verify_chirality.txt](results/verify_chirality.txt) | exact (Python integers) |
| Same polynomial identity on 715 interpolation points modulo four primes | [code/independent_tensor_check.py](code/independent_tensor_check.py) | [results/independent_tensor_check.txt](results/independent_tensor_check.txt) | exact modular arithmetic (NumPy int64) |
| 2×4 map from the HK 2016 formula at (2,2,2,1); tr([X,Y]^3) = -840i; real control 0 | [code/verify_qubit_ququart.py](code/verify_qubit_ququart.py) | [results/verify_qubit_ququart.txt](results/verify_qubit_ququart.txt) | exact |
| Block positivity of W and F (cited theorems) | [code/check_block_positivity.py](code/check_block_positivity.py) | [results/check_block_positivity.txt](results/check_block_positivity.txt) | floating point, seeded, tolerance -1e-9; sanity check only, not a proof |

Numbers quoted on the release page: 8.5×10^-10 (filtered-real lower bound), 8×10^-10 (ball radius), about 1.6×10^-7
(noise threshold η* = 1/6400000), -8/19199997 (witness value), dimensions 3×3 and 2×4, total dimension eight, degree
nine. All proved and checked exactly.

## Limits

The explicit state is extremely close to separable (ball 8×10^-10, noise threshold about 1.6×10^-7): a mathematical
statement, not an experimental proposal. W and D are Ha and Kye's. Whether part (ii) fails on 2×4 is open. Block
positivity rests on the cited Ha-Kye theorems (conditions checked exactly), with a floating-point cross-check only.
Experiments can still detect such states with several real measurements and classical post-processing. Open
problems: Section 9.

## Replay and records

Replay commands and comparison policy: [README.md](README.md), [ENVIRONMENT.txt](ENVIRONMENT.txt),
[REPLAY_RECEIPT.md](REPLAY_RECEIPT.md), CI [.github/workflows/replay.yml](.github/workflows/replay.yml) with
[.github/ci/compare_output.py](.github/ci/compare_output.py). Review: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Prior art
and limitation triage: [RESEARCH_GATES.json](RESEARCH_GATES.json). Bibliography: [refs.bib](refs.bib). Exact data:
[state.json](state.json). Licences: [LICENSE.md](LICENSE.md). Citation: [CITATION.cff](CITATION.cff). Hashes:
[MANIFEST.sha256](MANIFEST.sha256).
