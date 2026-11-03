exec(open('ctrl_z.py').read().split("from zstage import run")[0])
import xstage, sys
sys.path.insert(0, here)
from full import prefix
from sympl import setup
dmin2, syn2, Gam, b = setup(Cg, f)
basis, cost = prefix(m)
Zi = []
phi = {}
for x, z in Sp:
    phi.setdefault(prefcoords(x, z), syn2(tosuf(x & ~Tm, z & ~Tm)))
xs = [phi[1 << (m + i)] for i in range(m - 1)]
order = sorted(range(m - 1), key=lambda i: xs[i])
print('control own leaf (sites as given, unsorted is fine for the X-stage):',
      xstage.screen_leaf(xs, dmin2, Gam, b, m, cost, basis[:m-1], basis[m-1:]))
xstage.D = 5
print('control full run:', xstage.run(Cg, f, m))
