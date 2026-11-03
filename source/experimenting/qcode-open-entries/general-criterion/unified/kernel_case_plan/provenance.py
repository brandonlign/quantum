"""Check transitive saved census provenance without rerunning classifications."""
import json,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;sys.path.insert(0,str(ROOT))
from mixed_frames.stage3.common import atomic_json,sha

def main():
 checked={};reports=[]
 def pin(p):
  p=p.resolve();key=str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)
  checked[key]=sha(p)
 def verify(p,status):
  data=json.loads(p.read_text());assert data['status']==status,(p,data['status'])
  for name,digest in data['hashes'].items():
   path=Path(name);path=path if path.is_absolute() else ROOT/path
   assert sha(path)==digest,path;pin(path)
  pin(p);reports.append(str(p.relative_to(ROOT)));return data
 plan=verify(HERE/'plan.json','COMPLETE_PLAN_RELATIVE_TO_CLASSIFICATIONS')
 profiles=json.loads((ROOT/'kernel_capacity/profiles.json').read_text())
 assert profiles['status'].startswith('PASS') or profiles['status'].startswith('COMPLETE'),profiles['status']
 for name,digest in profiles['hashes'].items():
  p=ROOT/name;assert sha(p)==digest,p;pin(p)
 parent=verify(ROOT/'suffix_mass_check/verification8.json','PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS')
 child=verify(ROOT/'suffix_mass_check/verification9.json','PASS_RELATIVE_TO_CERTIFIED_RELEASE_PARENTS')
 assert child['parent_verification_sha256']==sha(ROOT/'suffix_mass_check/verification8.json')
 assert child['source_classes']==32865 and child['every_star_certificate_and_class_multiplicity_matches']
 assert parent['children']==26299 and parent['labelled_mass_lhs']==parent['labelled_mass_rhs']=='1593099450132480'
 for p in [Path(__file__),HERE/'PROOF.md',ROOT/'suffix_mass_check/PROOF.md',ROOT/'opus_workers/residual_plan_audit/AUDIT.md',ROOT/'opus_workers/suffix_mass_audit/AUDIT.md']:pin(p)
 assert plan['prefix_counts']=={'1,0':4,'2,0':6,'2,2':2,'3,0':7,'3,2':5}
 atomic_json(HERE/'provenance.json',dict(status='PASS_TRANSITIVE_SAVED_PROVENANCE',reports=reports,hashes=checked,scope='Verifies saved source/report digest closure; inherits explicitly documented certified parent completeness. No new classification or quantum extension exclusion.'))
 print(json.dumps(dict(status='PASS_TRANSITIVE_SAVED_PROVENANCE',pinned_files=len(checked))),flush=True)

if __name__=='__main__':main()
