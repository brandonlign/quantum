"""Soundness stress: m8 positives pinned to first 3 labels (many depth-4 nodes,
nonzero leaves). Frozen vs mitm4 byte-identical, bounded by --seconds 9."""
import json,os,sys,hashlib
from pathlib import Path
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
import fixed_kernel_mitm4.gate1 as G
rows=[];orig=G.compare
def cmp(name,model,data,c,extra,rs,**kw):
 if not name.endswith('prefix7'):return
 pins=extra[1].split(',')
 for k in (3,2):
  real=['--prefix',','.join(pins[:k]),'--seconds','9']
  o0,s0=G.run(G.FROZEN,model,data,c,real);o1,s1=G.run(G.NEW,model,data,c,real+['--mitm4',G.CAP,'--mitm4-check',G.CHECK])
  want='LEAF '+' '.join(pins)
  r=dict(name=name[:-1]+str(k),frozen_status=s0['status'],status=s1['status'],frozen_leaves=s0['leaves'],leaves=s1['leaves'],frozen_nodes=s0['nodes'],nodes=s1['nodes'],frozen_seconds=s0['seconds'],seconds=s1['seconds'],
   identical=(o0==o1) if s0['status']==s1['status']=='COMPLETE' else None,known_leaf_in_mitm=want in o1.splitlines(),**{k2:s1[k2] for k2 in ('mitm_tested','mitm_pruned','mitm_skipped','mitm_overflow','mitm_checked')})
  if s0['status']==s1['status']=='COMPLETE':assert o0==o1,r
  rs.append(r);print(json.dumps(r),flush=True)
  if s0['status']=='COMPLETE':break
G.compare=cmp;G.large(rows)
(HERE/os.environ.get('OUT','gate1_pin3.json')).write_text(json.dumps(dict(cap=int(G.CAP),check=int(G.CHECK),rows=rows),indent=1)+'\n')
