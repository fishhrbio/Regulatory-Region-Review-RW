"""Compare exact scientific outputs against the committed reference results."""
import hashlib
from pathlib import Path
import sys

reference=Path(__file__).resolve().parents[1]
rerun=Path(sys.argv[1]).resolve()
files=['results.json','qc.csv','state_coverage.csv','class_overlap.csv','paired_state_matrix.csv','paired_intervals.tsv.gz']
for name in files:
    a=(reference/'data/processed'/name).read_bytes()
    b=(rerun/'data/processed'/name).read_bytes()
    if hashlib.sha256(a).digest()!=hashlib.sha256(b).digest():
        raise AssertionError(f'Scientific output changed: {name}')
if (reference/'docs/index.html').read_bytes()!=(rerun/'docs/index.html').read_bytes():
    raise AssertionError('Generated explorer changed')
print('All six scientific data products and the generated explorer reproduce byte-for-byte.')
