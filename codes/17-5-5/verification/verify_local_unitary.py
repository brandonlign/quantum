"""Support-purity certificate against local-unitary equivalence to GF(4)-linearity."""
import json
from pathlib import Path
from fractions import Fraction
rows=json.loads((Path(__file__).parent/'generators.json').read_text())['rows']
words={0}
for row in rows:
    words |= {v^row for v in words}
if len(words)!=4096:
    raise RuntimeError('rank mismatch')
mask=5852
support=lambda v:(v&((1<<17)-1))|(v>>17)
shortened=sorted(v for v in words if support(v)&~mask==0)
if shortened != [0,221778000] or mask.bit_count()!=8:
    raise RuntimeError('shortening certificate mismatch')
purity=Fraction(len(shortened),1<<mask.bit_count())
print('PASS: subset', [i for i in range(17) if mask>>i&1])
print('Shortened stabilizer exactly',shortened)
print('Reduced normalized-code-projector purity',purity)
print('2^|T| times purity = 2, not a power of 4; impossible for GF(4)-linear stabilizers.')
