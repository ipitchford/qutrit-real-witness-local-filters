# Review summary (public, 0.1.0-candidate)

This records the internal review of the candidate: who looked, what they raised, and what was done. It is
AI-assisted internal review, not peer review or independent reproduction.

## What ran

1. **Panel (round 1).** Five reviewer roles, editor, methodology (proof audit), domain (literature, attribution and
   novelty), perspective (readability for non-specialists and significance) and devil's advocate, run as one combined
   Claude pass. The methodology reviewer re-derived the main chain by hand and with independent NumPy computations.
2. **Lead proof audit and clean replay** of the source bundle's programs from a fresh copy.
3. **One revision round** with extensions: a criterion in every dimension (Theorem 6), the qubit-ququart case promoted
   to a theorem (Theorem 10), a dimension-specificity remark, a rebuilt context-first introduction, definitions, a
   figure, a scope and open-problems section, and rewritten verification programs with explicit failure exits.
4. **Verification re-review and proof audit** of every new and rewritten proof (line by line, with the cited source
   statements compared and all five programs run under `python -O`): no critical or major issue; the minor wording,
   attribution and notation items were applied.
5. **Clean release replay** by the publisher from fresh copies, with and without `python -O`, plus a negative control.

No second-model (Astra/Codex) adversarial review took place for this release.

## Sources

Shen, Chen and Bian (Conjecture 9, Definition 1, Lemma 10, Theorems 8 and 11, Example 12, Section III), Ha and Kye 2012
(Theorem 2.2 and the indecomposability of the maps), Ha and Kye 2013 (Theorem 3.2, exposedness; the separable state;
the pairing convention), Ha and Kye 2016 (Theorem 3.2 and the map formula), Gurvits and Barnum (Theorem 1), Horodecki,
Horodecki and Horodecki (1996), Skowronek (2016), Wang et al. (Example 2), Fraser et al. (abstract), McKague, Mosca and
Gisin, Renou et al., Hickey and Gour, and Gühne and Tóth were obtained and checked against the statements the paper uses.
Woronowicz is credited only through Ha and Kye's own attribution. Straszewicz's theorem is not cited: the criterion has
a self-contained proof.

## Main objections and dispositions

| Objection | Disposition |
|---|---|
| Correctness | No mathematical error found by the panel, the lead's audit or the verification re-review. |
| Block positivity of the witness was asserted without a citation, and a Gram computation was described as an independent certificate | Ha and Kye (2012), Theorem 2.2, cited; its conditions are checked exactly in code (they hold with equality); the Gram computation is described as proving uniqueness given block positivity; a floating-point cross-check was added and labelled as not a proof. No exact sum-of-squares certificate was added (the cited theorem suffices). |
| No context for non-specialists | Introduction rebuilt: the problem, why it matters, the conjecture stated with SCB's sets E^LU and E^SLOCC, what was known, what is new, and a roadmap; abstract rewritten. |
| Definitions and notation | Definitions added (block positivity, witness, Choi matrix, decomposability, extreme and exposed points, filtered real); overloaded symbols renamed; every section opens with its purpose. |
| The programs did not check every certified constant and used bare `assert` | `verify_witness.py` rewritten to rebuild W and the separable state from Ha and Kye's data and check every constant of Section 6, the noise threshold and the state file; all checks exit with an error, also under `python -O` (negative control). The 2×4 program rebuilds the map from Ha and Kye's general formula and has a real control. |
| The two-copy statement was overframed | Now Proposition 11, framed as an ancilla resource: the second copy enters only through the local marginal and the filter is not unitary; one fresh ancilla qubit suffices with local unitaries; the copy number with local unitaries only is listed as open; contrasted with Wang et al.'s crossed-copy tests. |
| Attribution and literature | Ha and Kye 2012, McKague-Mosca-Gisin, Hickey-Gour and Gühne-Tóth added; Ha and Kye 2013 credited for the mechanism of Lemma 4; the witness and the separable state credited to Ha and Kye. |
| Significance and margins | An explicit noise-scale paragraph: ball radius 8×10^-10 and noise threshold about 1.6×10^-7, far too fragile for an experiment; open problems section. |
| The 2×4 case was undersold | Promoted to Theorem 10: part (i) fails in total dimension eight, the minimum. |

## Limits that remain

Whether part (ii) fails for a qubit and a ququart; for which dimensions an embedded two-qutrit counterexample stays
undetected; the largest noise tolerance among counterexamples; whether degree nine is the smallest separating degree;
the copy number with local unitaries only (Section 9 of the paper). Block positivity rests on the cited theorems with a
floating-point cross-check. The literature search is recorded in RESEARCH_GATES.json.
