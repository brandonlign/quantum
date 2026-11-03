"""Paired bounded fixed-first-image timings; no target exclusion claim."""
import json,os,signal,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from model import ExtensionModel
from fixed_suffix_search import quotient_input
from fixed_kernel.common import kernel_contexts,payload
from gauge_rankone_h5.common import normalize_word
from mixed_frames.stage3.common import atomic_json

def main():
 rows=json.loads((ROOT/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes'];C=rows[0]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C);data=quotient_input(model,[]);records=[]
 for pi,r in enumerate(json.loads((ROOT/'kernel_capacity/profiles.json').read_text())['records'][:4]):
  H=[normalize_word(h,r['p'],8) for h in r['H']];H,order,cs=kernel_contexts(H,8,data['Gamma'],with_group=True);c=cs[0];group=HERE/'benchmark_group.txt';group.write_text(f'{len(c["group"])} 7\n'+'\n'.join(' '.join(map(str,a)) for a in c['group'])+'\n')
  first=next(v for v,d in enumerate(data['distances']) if d>=5 and (v&c['radical_coeff'][0]).bit_count()%2==c['tags'][0]);base=[str(HERE/'z_enum_fixed'),'--gamma','--avoid-radical','--first-orbits','--kernel-data','--prefix-group',str(group),'--seconds','5']
  for xdom in [False,True]:
   start=time.monotonic();p=subprocess.run(base+(['--xdom'] if xdom else []),input=payload(model,data,[first],c),text=True,capture_output=True,check=True,timeout=20);s=json.loads(p.stderr);leaves=[list(map(int,l.split()[1:])) for l in p.stdout.splitlines()];assert len(leaves)==s['leaves']
   record=dict(profile=pi,first=first,xdom=xdom,stats=s,wall=time.monotonic()-start,leaves=leaves);records.append(record);print(json.dumps({k:v for k,v in record.items() if k!='leaves'}),flush=True)
 atomic_json(HERE/'xdom_benchmark.json',dict(status='BOUNDED_TIMINGS',records=records,scope='Only paired sample roots; projection/xdom not deployed in fixed-kernel production'))
if __name__=='__main__':main()
