"""Bounded real cross-producer audit and synthetic union-control checks."""
import copy
import json
import sys
from pathlib import Path
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import oracle_balance.aggregate as engine
from mixed_frames.stage3.common import atomic_json,sha


def main():
    closed='CONDITIONALLY_EXCLUDED_RESIDUAL_BRANCHES'
    real=engine.aggregate([0,32864])
    assert real['status']=='OPEN' and not real['errors'] and not real['orphans']
    assert len(real['records'])==2 and all(r['status']==closed for r in real['records'])
    template=dict(records=[dict(cls=0,status=closed),dict(cls=32864,status='OPEN')],
        errors=[],orphans=[],record_hashes={},counts={})
    vm=dict(template,records=[dict(cls=32864,status=closed)])
    controls=[]
    for mode in ['separate_complete','partial_only','mac_complete','duplicate_complete','vm_error','mac_error','orphan']:
        m,v=copy.deepcopy(template),copy.deepcopy(vm)
        if mode=='partial_only':v['records'][0]['status']='OPEN'
        if mode in ('mac_complete','duplicate_complete'):m['records'][1]['status']=closed
        if mode=='mac_complete':v['records'][0]['status']='OPEN'
        if mode=='vm_error':v['errors']=[dict(error='CONTROL_FOUND')]
        if mode=='mac_error':m['errors']=[dict(error='CONTROL_FOUND')]
        if mode=='orphan':v['orphans']=['CONTROL_ORPHAN']
        with patch.object(engine.original,'aggregate',side_effect=[m,v]),patch.object(engine,'atomic_json'):
            out=engine.aggregate([0,32864])
        assert len(out['records'])==2
        if mode in ('vm_error','mac_error','orphan'):assert out['status']=='ERROR'
        else:assert out['status']=='OPEN'  # A subset can never close the census.
        if mode=='partial_only':assert out['records'][1]['status']=='OPEN'
        if mode=='mac_complete':assert out['records'][1]['closure_producer']=='mac'
        if mode=='separate_complete':assert out['records'][1]['closure_producer']=='oracle'
        if mode=='duplicate_complete':assert out['counts'][closed]==2
        controls.append(dict(mode=mode,status='PASS'))
    atomic_json(ROOT/'oracle_balance/verification.json',dict(status='PASS_BOUNDED_UNION_ADMISSION',
        real_classes=[0,32864],real_report_sha256=sha(ROOT/'oracle_balance/aggregate.json'),
        controls=controls,hashes={str(p.relative_to(ROOT)):sha(p) for p in
            [Path(__file__),ROOT/'oracle_balance/aggregate.py',ROOT/'oracle_balance/watch.py',ROOT/'oracle_balance/stop_remote.py',ROOT/'oracle_balance/PROOF.md']},
        scope='Two real closed classes plus synthetic routing/failure controls. Not a complete target search or claim.'))
    print(json.dumps(dict(status='PASS_BOUNDED_UNION_ADMISSION',real_closed_classes=2,synthetic_controls=len(controls))),flush=True)


if __name__=='__main__':
    main()
