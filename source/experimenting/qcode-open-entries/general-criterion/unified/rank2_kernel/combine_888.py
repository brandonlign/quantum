"""Verified 888 coverage: new hull2/6 full-H sweep plus frozen hull0 proofs."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path[:0]=[str(ROOT),str(ROOT/'m8r1')]
from rank2_kernel.run import validate,fingerprints
from mixed_frames.stage4.run import read_manifest,verify_partition
from mixed_frames.stage4.common import HERE as OLD
from mixed_frames.stage3.common import sha,atomic_json

def main():
 plan,_=validate();new=json.loads((HERE/'aggregate.json').read_text())
 assert new['status']=='CONDITIONALLY_COMPLETE' and new['classes']==new['closed_classes']==41
 assert new['hashes']==fingerprints() and new['aggregation_source_sha256']==sha(HERE/'aggregate.py')
 assert len(new['record_hashes'])==1375
 for path,digest in new['record_hashes'].items():assert sha(ROOT/path)==digest,path
 manifest=read_manifest();h0={r['cls'] for r in plan['records'] if r['q']==0};assert len(h0)==15
 old_hashes={};cases=[];missing=[]
 for row in manifest['classes']:
  if row['class_index'] not in h0:continue
  total=0
  for first in row['first_images']:
   folder=OLD/'runs'/row['key']/str(first);path=folder/'checkpoint.json'
   if not path.exists() or json.loads(path.read_text())['status']!='COMPLETE':missing.append(dict(key=row['key'],first=first));continue
   cp=verify_partition(folder,row,first,manifest['hashes']);assert cp['excluded'] and cp['surviving_solutions']==0
   old_hashes[str(path.relative_to(ROOT))]=sha(path);total+=1
  cases.append(dict(key=row['key'],complete=total==len(row['first_images']),verified_roots=total,roots=len(row['first_images'])))
 assert len(cases)==65
 hashes={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'aggregate.py',HERE/'aggregate.json',HERE/'REVIEW.md',OLD/'manifest.json',OLD/'exclusions.json',OLD/'run.py',OLD/'PROOF.md']}
 hashes.update(manifest['hashes'])
 result=dict(status='CONDITIONALLY_COMPLETE_888_COLUMN' if not missing else 'OPEN',hull2_hull6_classes=41,hull0_classes=15,classes=56,new_slices=1375,new_Y_leaves_rejected=new['Y_leaves_completely_rejected'],legacy_hull0_verified_roots=len(old_hashes),legacy_cases=cases,missing=missing,hashes=hashes,legacy_checkpoint_hashes=old_hashes,new_record_hashes=new['record_hashes'],scope='Complete exact rank-two 888 column conditional on recorded classifications and prior recentering/purity/profile certificates; rank-one remains open, not a whole-target exclusion')
 atomic_json(HERE/'combined_888.json',result)
 print(json.dumps({k:v for k,v in result.items() if k not in ('legacy_cases','missing','hashes','legacy_checkpoint_hashes','new_record_hashes')}),flush=True)
if __name__=='__main__':main()
