"""Exact equivalence grouping of complete length9 extension representatives."""
from pathlib import Path
from collections import Counter
import json,hashlib
from pynauty import Graph,certificate
BASE=Path(__file__).resolve().parent

def cert(gens,n):
    words=[0]
    for g in gens:words += [w^g for w in words]
    count=len(words);adj={i:set() for i in range(count+3*n)}
    for row,w in enumerate(words):
        for j in range(n):
            s=((w>>j)&1)|(((w>>(n+j))&1)<<1)
            if s:adj[row].add(count+3*j+s-1)
    for j in range(n):
        for s in range(3):
            adj[count+3*j+s].update(count+3*j+t for t in range(3) if t!=s)
    return certificate(Graph(count+3*n,adjacency_dict=adj,
                       vertex_coloring=[set(range(count)),set(range(count,count+3*n))]))

def main():
    data=BASE/'suffix9_extensions.json';source=json.loads(data.read_text());classes={}
    for i,item in enumerate(source['lifts']):
        key=cert(item['generators'],9).hex()
        if key not in classes:classes[key]=dict(representative=item,extension_indices=[])
        classes[key]['extension_indices'].append(i)
    hist=Counter(v['representative']['hull'] for v in classes.values())
    vals=[]
    for key,item in sorted(classes.items()):
        vals.append(dict(class_index=len(vals),certificate_hex=key,**item))
    report=dict(status='COMPLETE',equivalence='site permutation and GL(2,2) per site',
        extension_representatives=len(source['lifts']),inequivalent_classes=len(vals),
        hull_counts=dict(hist),remaining_target_classes=sum(v['representative']['hull']==0 for v in vals),
        classes=vals,scope='Exact graph equivalence quotient of complete parent-extension enumeration; '
                          'not a full stabilizer extension search.',
        sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),data]})
    (BASE/'suffix9_classes.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['classes','sha256']}))
if __name__=='__main__':main()
