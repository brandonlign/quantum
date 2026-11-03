"""Cache only the immutable local frame in a separate Y solver module."""
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
s=(ROOT/'kernel_y_parity/solve.py').read_text()
old='from fixed_kernel.common import standard,kernel'
new='from fixed_kernel.common import standard as build_frame,kernel\nfrom functools import lru_cache\nstandard=lru_cache(maxsize=8)(build_frame)'
assert s.count(old)==1
(HERE/'solve.py').write_text(s.replace(old,new))
