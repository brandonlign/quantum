"""Read-only scan of mixed_frames/stage4 checkpoints (no verification, no writes there)."""
import json,sys,statistics
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parent.parent;S4=ROOT/'mixed_frames/stage4'
man=json.loads((S4/'manifest.json').read_text())
st=defaultdict(int);byprof=defaultdict(lambda:defaultdict(int));secs=defaultdict(list);leaves=0;surv=0;nodes=[];running=[]
for row in man['classes']:
    prof=''.join(map(str,row['profile']))
    for first in row['first_images']:
        p=S4/'runs'/row['key']/str(first)/'checkpoint.json'
        if not p.exists():s='PENDING'
        else:
            try:c=json.loads(p.read_text())
            except Exception:s='UNREADABLE'
            else:
                s=c['status']
                if s=='COMPLETE':
                    secs[prof].append(c['seconds']);leaves+=c['leaves'];surv+=c['surviving_solutions']
                    nodes.append((c['enumerator']['nodes'],row['key'],first,c['seconds']))
                if s=='RUNNING':running.append((row['key'],first,c['seconds'],c['leaves']))
        st[s]+=1;byprof[prof][s]+=1
out=dict(total=man['partitions'],status=dict(st),by_profile={k:dict(v) for k,v in byprof.items()},
  complete_leaves=leaves,surviving_solutions=surv,running=running,
  complete_seconds={k:dict(n=len(v),sum=round(sum(v),1),median=round(statistics.median(v),4),max=round(max(v),1)) for k,v in secs.items()},
  top_nodes=sorted(nodes,reverse=True)[:10],total_nodes=sum(n for n,*_ in nodes))
print(json.dumps(out,indent=1))

def per_class():
    rows=defaultdict(lambda:dict(done=0,pending=0,sec=0.0,nodes=0,leaves=0))
    same=diff=0;cmp=defaultdict(dict)
    for row in man['classes']:
        prof=''.join(map(str,row['profile']));k=(row['class_index'],prof)
        for first in row['first_images']:
            p=S4/'runs'/row['key']/str(first)/'checkpoint.json'
            c=json.loads(p.read_text()) if p.exists() else None
            if c and c['status']=='COMPLETE':
                r=rows[k];r['done']+=1;r['sec']+=c['seconds'];r['nodes']+=c['enumerator']['nodes'];r['leaves']+=c['leaves']
                if prof in('10','11','20'):cmp[(row['class_index'],first)][prof]=(c['enumerator']['nodes'],c['leaves'])
            else:rows[k]['pending']+=1
    for v in cmp.values():
        if len(v)==3:
            if len(set(v.values()))==1:same+=1
            else:diff+=1
    return rows,same,diff
if __name__=='__main__' and 'classes' in sys.argv:
    rows,same,diff=per_class()
    print('profile 10/11/20 identical (nodes,leaves):',same,'different:',diff)
    for (i,pf),r in sorted(rows.items()):print(i,pf,r)
