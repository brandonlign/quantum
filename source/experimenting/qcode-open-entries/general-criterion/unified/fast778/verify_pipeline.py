"""Actual streaming worker and independent saved-evidence rejection tests."""
import copy
import json
from pathlib import Path
import tempfile
from model import ExtensionModel
from mixed_frames.frames import adapted_frame
from mixed_frames.stage3.common import input_text
from . import run
from .common import HERE, hashes, sha, atomic_json


def main():
    m=ExtensionModel(6,0,3,4,[15,240],[])
    f=adapted_frame(4,m.R)
    original_here, original_target=run.HERE,run.target
    expected=hashes()
    try:
        with tempfile.TemporaryDirectory(dir=HERE,prefix='pipeline_') as name:
            dest=Path(name)
            run.HERE=dest
            run.target=lambda frame,index:(m,f)
            path=dest/'input.txt'
            path.write_text(input_text(m,f,3))
            row=dict(key='fixture',frame='fixture',class_index=0,profile=[0,0],budget=None,
                     d=3,mu_out=3,input='input.txt',input_sha256=sha(path))
            result=run.partition((row,3,expected))
            assert result['status']=='COMPLETE' and result['surviving_solutions']>0 and not result['excluded'],result
            folder=dest/'runs/fixture/3'
            assert run.verify_partition(folder,row,3,expected)==result
            corrupt=copy.deepcopy(result)
            corrupt['surviving_solutions']=0;corrupt['excluded']=True
            atomic_json(folder/'checkpoint.json',corrupt)
            try:run.verify_partition(folder,row,3,expected)
            except AssertionError:pass
            else:raise AssertionError('False zero-survivor evidence accepted')
            atomic_json(folder/'checkpoint.json',result)
            (folder/'enumerator.jsonl').write_text('')
            try:run.verify_partition(folder,row,3,expected)
            except AssertionError:pass
            else:raise AssertionError('Missing terminal evidence accepted')
            report=dict(status='PASS',hashes=expected,leaves=result['leaves'],
                        surviving_solutions=result['surviving_solutions'],
                        false_exclusion_rejected=True,missing_summary_rejected=True)
    finally:run.HERE,run.target=original_here,original_target
    atomic_json(HERE/'pipeline_verification.json',report)
    print(json.dumps({k:v for k,v in report.items() if k!='hashes'}))


if __name__=='__main__':main()
