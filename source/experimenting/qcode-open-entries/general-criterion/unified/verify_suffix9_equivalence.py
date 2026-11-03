"""Different graph encoding for exact suffix equivalence partition replay."""
from pathlib import Path
import json,hashlib
from pynauty import Graph,certificate
BASE=Path(__file__).resolve().parent

def star_certificate(gens,n):
    words=[0]
    for g in gens:words += [w^g for w in words]
    W=len(words);S=W;G=W+3*n;V=G+n;adj={i:set() for i in range(V)}
    def connect(a,b):adj[a].add(b);adj[b].add(a)
    for row,w in enumerate(words):
        for j in range(n):
            s=((w>>j)&1)|(((w>>(n+j))&1)<<1)
            if s:connect(row,S+3*j+s-1)
    for j in range(n):
        for s in range(3):connect(G+j,S+3*j+s)
    graph=Graph(V,adjacency_dict=adj,vertex_coloring=[set(range(W)),set(range(S,G)),set(range(G,V))])
    return certificate(graph)

def main():
    raw=json.loads((BASE/'suffix9_extensions.json').read_text())['lifts']
    saved=json.loads((BASE/'suffix9_classes.json').read_text())['classes'];groups={}
    for i,v in enumerate(raw):groups.setdefault(star_certificate(v['generators'],9),[]).append(i)
    independent={tuple(v) for v in groups.values()};original={tuple(v['extension_indices']) for v in saved}
    assert independent==original
    report=dict(status='PASS',extension_representatives=len(raw),inequivalent_classes=len(groups),
                identical_partitions=True,independent_encoding='codeword/symbol/coordinate vertices in three colors, stars',
                original_encoding='codeword/symbol vertices in two colors, coordinate triangles',
                scope='Independent graph encoding; same pynauty canonicalization library and parent census.',
                sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE/'suffix9_classes.json',BASE/'suffix9_extensions.json']})
    (BASE/'suffix9_equivalence_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'}))
if __name__=='__main__':main()
