"""Check syndrome quotient conversion independently on actual target data."""
import json
import random
from pathlib import Path
from run_hull5_planes import quotient_input,heavy_planes,ExtensionModel,z_condition
from model import gram,lift_coordinate,reduced_basis,symp,span,weight

HERE=Path(__file__).resolve().parent
classes=json.loads((HERE.parent/'opus_parallel/suffix9_dimension5/suffix9_dim5.json').read_text())['classes']
rng=random.Random(1736);records=[]
for index in [150,181,217]:
    C=classes[index]['representative']['generators'];model=ExtensionModel(17,3,6,8,[255<<8],C)
    planes=heavy_planes(model.C);assert len(planes)==1
    data=quotient_input(model,planes[0]);labels=data['labels']
    assert len(set(labels))==2048
    expected={s for s in range(8192) if all((s&n).bit_count()%2==0 for n in data['Ncoords'])}
    assert set(labels)==expected
    for small,old in enumerate(labels):
        assert data['distances'][small]==min(weight(model.representatives[old]^c,9) for c in span(model.C))
    # Compare direct restricted-Gram compatibility in both coordinate systems.
    for _ in range(200):
        small=reduced_basis([rng.randrange(2048) for _ in range(7)])
        if len(small)!=7:continue
        old=[lift_coordinate(s,data['units']) for s in small]
        assert z_condition(small,data['Gamma'],11)==z_condition(old,model.Gamma,13)
    records.append(dict(cls=index,planes=len(planes),labels_checked=2048,
                        physical_coset_words=2048*32,random_compatibility_trials=200))
result=dict(status='PASS',records=records,
            scope='Exact annihilator bijection, physical coset minima, and original-vs-quotient compatibility')
(HERE/'heavy_plane_quotient_verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
