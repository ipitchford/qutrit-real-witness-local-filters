# SPDX-License-Identifier: MIT
"""Compare fresh outputs with the recorded outputs in results/ under the declared policy.

usage: compare_output.py NEW_DIR RECORDED_RESULTS_DIR

NEW_DIR holds <script>.txt for the five programs and verification_witness.json written by the fresh run.
Policy (declared before tagging; the same in README.md, ENVIRONMENT.txt and REPLAY_RECEIPT.md):
- verify_chirality, verify_qubit_ququart, independent_tensor_check: identical line by line.
- verify_witness: identical line by line except the float display lines "alpha:" and
  "lower bound on filtered real expectations:", whose values agree within 1e-9 relative.
- verification_witness.json: every field identical except the float display fields tau_r, cubic_change_bound, alpha,
  filtered_real_expectation_lower_bound, psd_margin and ppt_margin, which agree within 1e-9 relative.
- check_block_positivity (floating-point sanity check): same lines with the numbers masked; each printed minimum
  >= -1e-9; final line identical.
- Interpreter version fields are not part of any output and are not compared.
Exits 1 on any difference.
"""
import json
import re
import sys
from pathlib import Path

REL = 1e-9
EXACT = ['verify_chirality', 'verify_qubit_ququart', 'independent_tensor_check']
FLOAT_LINES = ('alpha:', 'lower bound on filtered real expectations:')
FLOAT_FIELDS = {'tau_r', 'cubic_change_bound', 'alpha', 'filtered_real_expectation_lower_bound', 'psd_margin', 'ppt_margin'}
NUM = re.compile(r'-?\d+(?:\.\d+)?(?:e[-+]?\d+)?')


def close(a, b):
    return abs(a - b) <= REL * max(abs(a), abs(b))


def lines(p):
    return Path(p).read_text(encoding='utf-8').splitlines()


def main():
    new, old = Path(sys.argv[1]), Path(sys.argv[2])
    bad = []
    for s in EXACT:
        if lines(new / f'{s}.txt') != lines(old / f'{s}.txt'):
            bad.append(f'{s}: output differs')
    a, b = lines(new / 'verify_witness.txt'), lines(old / 'verify_witness.txt')
    if len(a) != len(b):
        bad.append('verify_witness: line count differs')
    for x, y in zip(a, b):
        if x == y:
            continue
        if x.startswith(FLOAT_LINES) and x.split(':')[0] == y.split(':')[0]:
            if not close(float(x.split(':')[1]), float(y.split(':')[1])):
                bad.append(f'verify_witness: {x!r} vs {y!r} beyond 1e-9 relative')
        else:
            bad.append(f'verify_witness: {x!r} != {y!r}')
    if not a or a[-1] != 'ALL 120 CHECKS PASSED':
        bad.append('verify_witness: final line')
    ja, jb = (json.loads(Path(p).read_text()) for p in (new / 'verification_witness.json', old / 'verification_witness.json'))
    if set(ja) != set(jb):
        bad.append('json: key sets differ')
    for k in sorted(set(ja) & set(jb)):
        if k in FLOAT_FIELDS:
            if not (isinstance(ja[k], float) and close(ja[k], jb[k])):
                bad.append(f'json {k}: {ja[k]} vs {jb[k]}')
        elif ja[k] != jb[k]:
            bad.append(f'json {k} differs')
    a, b = lines(new / 'check_block_positivity.txt'), lines(old / 'check_block_positivity.txt')
    if [NUM.sub('#', x) for x in a] != [NUM.sub('#', x) for x in b]:
        bad.append('check_block_positivity: lines differ outside numbers')
    for x in a[:-1]:
        v = float(x.rsplit(':', 1)[1])
        if v < -1e-9:
            bad.append(f'check_block_positivity: {x!r} below -1e-9')
    if not a or a[-1] != b[-1]:
        bad.append('check_block_positivity: final line')
    if bad:
        print('DIFFERENCES:')
        for m in bad:
            print(' ', m)
        sys.exit(1)
    print('OK: exact outputs identical; verify_witness and JSON identical outside float display fields (within 1e-9 '
          'relative); block-positivity minima >= -1e-9')


if __name__ == '__main__':
    main()
