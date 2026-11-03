"""Recompute the actual q=3 branch screens on every physical kernel plane.

Counts avoided fixed-kernel search tasks, not unmeasured runtime savings.
Classification representatives and orbit partitions are the frozen inputs.
Every plane is regenerated from the radical and checked against that partition.
"""
from pathlib import Path
import collections,hashlib,itertools,json,sys,time
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'experimenting/qcode-open-entries/general-criterion/unified'
sys.path.insert(0,str(BASE))
from model import span,gram,nullspace,lift_coordinate,weight

def screen(plane):
    ws=[((w&511)|(w>>9)).bit_count() for w in plane]
    support=[(w&511)|(w>>9) for w in plane]
    if sum(w<8 for w in ws)>2:return 'three_light_words'
    if any(ws[i]==8 and any(j!=i and support[j]&~support[i]==0 for j in range(3)) for i in range(3)):
        return 'rank_two_support_already_handled'
    return None

def main():
    started=time.monotonic();report={'status':'RUNNING','scope':__doc__,'branches':[]}
    sources=[('q3_t0_h5','rankone_recenter/plane_plan.json','representative',1,118),
             ('q3_t2_h3','fixed_kernel/hull3_plan.json','N',5,5995)]
    for name,path,basis_key,types,expected in sources:
        source=BASE/path;plan=json.loads(source.read_text());assert len(plan['records'])==expected
        counts=collections.Counter();physical=collections.Counter()
        for row in plan['records']:
            C=row['C']
            rad=[lift_coordinate(x,C) for x in nullspace(gram(C,9),len(C))]
            words=[w for w in span(rad) if w]
            allplanes={tuple(sorted((a,b,a^b))) for a,b in itertools.combinations(words,2)}
            supplied=set()
            for orbit in row['orbits']:
                planes=[tuple(sorted(p)) for p in orbit['members']]
                assert len(planes)==len(set(planes))
                assert all(p in allplanes for p in planes)
                assert not supplied.intersection(planes);supplied.update(planes)
                reasons={screen(p) for p in planes};assert len(reasons)==1
                reason=reasons.pop() or 'remaining'
                assert tuple(sorted(span(orbit[basis_key])[1:])) in planes
                counts[reason]+=1;physical[reason]+=len(planes)
            assert supplied==allplanes,(name,row['cls'])
        before=sum(counts.values())*types;after=counts['remaining']*types
        result={'branch':name,'shortened_code_classes':expected,'local_types_per_unscreened_orbit':types,
                'physical_plane_counts':dict(physical),'orbit_counts':dict(counts),
                'search_tasks_before_screens_after_symmetry':before,'search_tasks_after_screens':after,
                'search_tasks_avoided':before-after,'percent_search_tasks_avoided':100*(before-after)/before,
                'input_file':path,'input_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
                'every_physical_plane_regenerated_and_screen_recomputed':True}
        report['branches'].append(result);print(json.dumps(result),flush=True)
    report.update(status='PASS',elapsed_seconds=time.monotonic()-started,
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (ROOT/'benchmarks/results/actual_17_screens.json').write_text(json.dumps(report,indent=2))
if __name__=='__main__':main()
