"""m8 rows only (large positives + target roots) at a chosen CAP/CHECK; writes OUT."""
import json,os,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import fixed_kernel_mitm4.gate1 as G
rows=[];G.large(rows)
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,identification_orbits
from gauge_rankone_h5.common import normalize_word
classes=json.loads((G.ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];profiles=json.loads((G.ROOT/'kernel_capacity/profiles.json').read_text())['records']
C=classes[0]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[])
for ri in range(4):
 row=profiles[ri];H,_,cs=kernel_contexts([normalize_word(h,row['p'],8) for h in row['H']],8,data['Gamma'],with_group=True);c=cs[identification_orbits(H,8,cs)['orbits'][0]['representative']]
 reps=[v for v,d in enumerate(data['distances']) if d>=5 and all(((v&r).bit_count()%2)==((c['tags'][0]>>j)&1) for j,r in enumerate(c['radical_coeff']))]
 G.compare(f'target_h1_c0_prof{ri}_first115',model,data,c,['--prefix','115'],rows,reps=reps,group='only')
(HERE/os.environ['OUT']).write_text(json.dumps(dict(cap=int(G.CAP),check=int(G.CHECK),rows=rows),indent=1)+'\n')
