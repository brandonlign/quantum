"""Strict residual q1/q2 aggregation, reconstructing every covered slice."""
import argparse,json,re,sys
from collections import Counter,defaultdict
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rankone_residual.common import validate,fingerprints,slice_setup,CHUNK,digest
from rankone_residual.search import expected_slices
from fixed_kernel.symmetry import suffix_slice
from model import reduced_basis
from mixed_frames.stage3.common import atomic_json,sha

def verify_chunk(path,row,s,k,hashes):
 r=json.loads(path.read_text());branch=s['branch'];q,t=s['q'],s['t'];oi=s['oi'];prefix=s['prefix'];N=branch['N_orbits'][oi]['N']
 assert r['hashes']==hashes,'semantic fingerprints'
 assert r['backend'] in ('compact','joint'),'unregistered native backend'
 if 'pool_hashes' in r:
  from rankone_residual.pool import pool_hashes
  assert r['pool_hashes']==pool_hashes(),'pool fingerprints'
 assert (r['cls'],r['hull'],r['q'],r['t'],r['orbit_index'],r['profile_index'],r['context_index'],r['chunk'])==(row['cls'],row['hull'],q,t,oi,prefix['profile_index'],s['ci'],k),'task identity'
 assert r['N']==N and r['N_members']==branch['N_orbits'][oi]['members'],'N coverage'
 model,data,H,order,contexts,io=slice_setup(tuple(row['C']),tuple(N),tuple(prefix['H']),prefix['p'],q,t)
 assert r['H']==H and r['site_order']==order,'normalized H'
 assert r['identifications']==len(contexts) and r['identification_orbit']==s['members'],'identification coverage'
 assert sorted(i for o in io['orbits'] for i in o['members'])==list(range(len(contexts)))
 ss=suffix_slice(model,N,data,contexts[s['ci']]);reps=ss['first_representatives'];mine=reps[k*CHUNK:(k+1)*CHUNK]
 assert r['orbit_digest']==s['orbit_digest']==digest(ss),'physical first-label orbits'
 assert r['chunks']==s['chunks']==max(1,-(-len(reps)//CHUNK)) and r['chunk_size']==CHUNK
 assert r['first']==mine and r['first_representatives_total']==len(reps),'first-label partition'
 if r['status']!='EXCLUDED_SLICE_CHUNK':
  assert r['status'] in ('LIMIT','FOUND'),r['status'];return r['status']
 assert r['z']['status']=='COMPLETE' and r['z'].get('reason') is None,'unfinished Z search'
 leaves=r['leaves'];assert r['z']['leaves']==len(leaves),'missing Y leaves'
 actual=Counter(e['result']['status']+':'+str(e['result']['reason']) for e in leaves)
 assert dict(actual)==r['tally'],'Y tally'
 xs=[tuple(e['X']) for e in leaves];assert len(xs)==len(set(xs)),'duplicate Z leaves'
 for e in leaves:
  x=e['X'];a=e['result'];assert len(x)==7 and x[0] in mine and len(reduced_basis(x))==7
  assert all(isinstance(v,int) and 0<=v<len(data['labels']) for v in x)
  assert a['status']=='COMPLETE' and a['complete'] and a['fixed_z_excluded'] and not a['solutions'],'unfinished or positive Y leaf'
  assert a['reason']!='search_limit','Y limit mislabeled'
 assert 'physical' not in r and 'candidate_file' not in r,'candidate cannot be an exclusion'
 return 'EXCLUDED_SLICE_CHUNK'

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--classes',type=int,nargs='+');ap.add_argument('--folder',type=Path,default=HERE/'runs');ap.add_argument('--output',type=Path,default=HERE/'aggregate.json');args=ap.parse_args();assert __debug__
 plan=validate(require_controls=True);hashes=fingerprints();lookup={r['cls']:r for r in plan['records']};selected=set(lookup) if args.classes is None else set(args.classes);assert selected<=set(lookup)
 present=defaultdict(set);orphans=[]
 for path in args.folder.glob('*.json'):
  m=re.fullmatch(r'(\d+)_q[12]t[02]_\d+_\d+_\d+_\d+\.json',path.name)
  if not m:orphans.append(str(path));continue
  cls=int(m[1])
  if cls not in lookup:orphans.append(str(path));continue
  if cls in selected:present[cls].add(path.name)
 records=[];counts=Counter();record_hashes={};errors=[]
 for cls in sorted(selected):
  row=lookup[cls];branch_counts=defaultdict(Counter);seen=set()
  if not present[cls]:
   records.append(dict(cls=cls,status='PENDING',branches=[dict(q=b['q'],t=b['t'],status='PENDING') for b in row['branches']]));counts['PENDING_CLASSES']+=1;continue
  for s in expected_slices(plan,row):
   q,t=s['q'],s['t'];ct=Counter()
   for k,name in enumerate(s['names']):
    seen.add(name);path=args.folder/name
    if not path.exists():ct['PENDING']+=1;continue
    try:
     status=verify_chunk(path,row,s,k,hashes);ct[status]+=1;record_hashes[str(path.resolve())]=sha(path)
    except Exception as e:ct['ERROR']+=1;errors.append(dict(cls=cls,path=str(path),error=repr(e)))
   counts.update(ct)
   closed=ct==Counter(EXCLUDED_SLICE_CHUNK=s['chunks'])
   branch_counts[q,t]['CLOSED_SLICES' if closed else 'OPEN_SLICES']+=1
  if present[cls]-seen:orphans.extend(str(args.folder/name) for name in sorted(present[cls]-seen))
  branches=[dict(q=b['q'],t=b['t'],status='CONDITIONALLY_EXCLUDED_BRANCH' if branch_counts[b['q'],b['t']]['CLOSED_SLICES']>0 and branch_counts[b['q'],b['t']]['OPEN_SLICES']==0 else 'OPEN',counts=dict(branch_counts[b['q'],b['t']])) for b in row['branches']]
  closed=all(b['status']=='CONDITIONALLY_EXCLUDED_BRANCH' for b in branches);counts['CLOSED_CLASSES' if closed else 'OPEN_CLASSES']+=1
  records.append(dict(cls=cls,status='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' if closed else 'OPEN',branches=branches))
 complete=selected==set(lookup) and len(records)==len(lookup) and all(r['status']=='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES' for r in records) and not errors and not orphans
 out=dict(status='CONDITIONALLY_COMPLETE_RESIDUAL_RANKONE' if complete else ('ERROR' if errors or orphans else 'OPEN'),classes=len(lookup),selected_classes=len(selected),counts=dict(counts),records=records,errors=errors,orphans=orphans,record_hashes=record_hashes,hashes=hashes,controls_sha256=sha(HERE/'controls.json'),aggregation_source_sha256=sha(Path(__file__)),scope='Only q1/t0 hull1, q2/t2 hull1 and q2/t0 hull3, conditional on recorded classifications and exact physical symmetry coverage. The separate q3 branches and whole target require their own proofs.')
 atomic_json(args.output,out);print(json.dumps({k:v for k,v in out.items() if k not in ('records','hashes','record_hashes')}),flush=True)

if __name__=='__main__':main()
