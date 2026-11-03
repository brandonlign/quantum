"""Exact inverse for any pure gauge base; original local rank is one."""
import json,subprocess
from pathlib import Path
from model import reduced_basis,span,weight,symp,dual,nullspace,lift_coordinate,embed
from gauge_fix import supported_subspace
HERE=Path(__file__).resolve().parent

def ordered_basis(S,n,m):
 S=reduced_basis(S);H=supported_subspace(S,n,range(m));p=((1<<m)-1)<<n
 B=[p]
 assert p in span(H)
 for w in H:
  if len(reduced_basis(B+[w]))>len(B):B.append(w)
 h=len(B)
 for w in S:
  if len(reduced_basis(B+[w]))>len(B):B.append(w)
 assert len(B)==len(S)
 return B,h

def select(S,n,m,d,cut,budget,seconds=None,nodes=None):
 B,h=ordered_basis(S,n,m);r=len(B)
 text=f'{n} {d} {m} {h} {cut} {budget} {r}\n'+' '.join(map(str,B))+'\n'
 cmd=[str(HERE/'inverse')]
 if seconds is not None:cmd+=['--seconds',str(seconds)]
 if nodes is not None:cmd+=['--nodes',str(nodes)]
 p=subprocess.run(cmd,input=text,capture_output=True,text=True,check=True);result=json.loads(p.stdout)
 assert result['status'] in ('FOUND','EXCLUDED_BASE','LIMIT')
 if result['status']=='FOUND':
  rows=result['rows'];assert len(rows)==h-1
  assert all((v&((1<<h)-1))==1<<(i+1) for i,v in enumerate(rows))
  out=reduced_basis([lift_coordinate(c,B) for c in nullspace(rows,r)])
  assert len(out)==r-h+1 and all(symp(a,b,n)==0 for a in out for b in out)
  D=dual(out,n);assert min(weight(w,n) for w in span(D) if w)>=d
  assert sum(weight(w,n)<cut for w in span(out) if w)<=budget
  assert supported_subspace(out,n,range(m))==[embed(((1<<m)-1)<<m,0,m,n-m)]
  result.update(generators=out,normalizer=D)
 return result
