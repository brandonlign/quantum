"""Continue the live fast-frame queue while reserving one rank-one core."""
import json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from scheduler.handoff_fast_frames import log,python_manager,stop_tree,launch

def main():
 log('waiting for live fast778; rank-one owns the paused fast688 core')
 while python_manager('fast778.run'):time.sleep(30)
 aggregate=json.loads((ROOT/'fast778/aggregate.json').read_text())
 assert aggregate['status']=='COMPLETE' and aggregate['verified']
 from fast778.common import validate_manifest as check778
 from speedup_review.common import validate_manifest as check688
 check778();check688()
 for pid in python_manager('speedup_review.run'):stop_tree(pid)
 lease=ROOT/'rankone_light_budget/core_lease.json'
 data=json.loads(lease.read_text());data['status']='RELEASED_OLD688_TREE_REPLACED';lease.write_text(json.dumps(data,indent=2)+'\n')
 p=launch('speedup_review.run',6,'speedup_review/run_688.log')
 log('6 fast688 + 1 existing 888 + 1 rank-one workers = 8 cores')
 code=p.wait()
 if code!=0:raise RuntimeError(f'fast688 exited {code}')
 log('fast688 finished; strict aggregate and independent checks still required')
if __name__=='__main__':main()
