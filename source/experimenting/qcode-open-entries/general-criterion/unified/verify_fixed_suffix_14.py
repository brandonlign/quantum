"""Reuse the generic quotient search on the localized impure fourteen-site reduction."""
import hashlib
import json
from pathlib import Path
from model import ExtensionModel,nullspace,lift_coordinate
from fixed_suffix_search import search


def main():
    here=Path(__file__).resolve().parent
    source=Path(__file__).resolve().parents[2]/'release16/census/cells.json'
    codes=json.loads(source.read_text())['10,10'];records=[]
    # Binary kernel distance>=6 on four sites forces q0, hence hull4
    # and the whole four-dimensional radical is the suffix stabilizer.
    for i,C in enumerate(codes):
        model=ExtensionModel(14,3,5,4,[15<<4],C)
        if model.hull!=4:
            r=dict(class_index=i,status='EXCLUDED_KERNEL_HULL',hull=model.hull)
        else:
            N=[lift_coordinate(v,model.P) for v in nullspace(model.Gamma,model.b)]
            assert len(N)==4
            r=search(model,N,4,seconds=5);r.update(class_index=i,hull=4)
        records.append(r)
    assert len(records)==37 and all(r['status'] in ['EXCLUDED_KERNEL_HULL','EXCLUDED_FIXED_N'] for r in records)
    result=dict(status='COMPLETE_LOCALIZED_14_REPLAY',records=records,
                hull_exclusions=sum(r['status']=='EXCLUDED_KERNEL_HULL' for r in records),
                generic_search_exclusions=sum(r['status']=='EXCLUDED_FIXED_N' for r in records),
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in
                        [source,Path(__file__),here/'fixed_suffix_search.py',here/'model.py',
                         here/'signed_extension.py',here/'m8r1/z_enum_avoid_radical.c',
                         here/'m8r1/z_enum_avoid_radical']},
                scope='Same generic fixed-N quotient and C runner on localized impure[[14,3,5]] reduction; shared existing census and search core; not a new independent paper proof')
    (here/'fixed_suffix_14_replay.json').write_text(json.dumps(result,indent=2)+'\n')
    print({k:v for k,v in result.items() if k not in ['records','sha256']})


if __name__=='__main__':main()
