"""Physical low-suffix words and necessary rank-one light-word counts."""
from pathlib import Path
from model import lift_coordinate, span, weight

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
BINARY=HERE/'z_enum_light'
OLD=ROOT/'m8r1/z_enum_symmetry'


def physical_data(model,N,data,cut):
    # A quotient coefficient is the evaluation of the original coefficient
    # vector on each of the quotient syndrome units. Multiple physical
    # suffix words with this coefficient differ by N and must each be counted.
    units=data['units'];entries=[];physical=[]
    for c in range(1<<model.b):
        w=lift_coordinate(c,model.P);wt=weight(w,model.f)
        if wt+model.m//2>=cut:continue
        a=sum(((c&u).bit_count()%2)<<j for j,u in enumerate(units))
        if a:
            entries.append((a,wt));physical.append(w)
    assert len(physical)==len(set(physical))
    baseline=sum(weight(r,model.m)+weight(n,model.f)<cut
                 for r in span(model.R) for n in span(N) if r or n)
    return dict(entries=entries,physical_words=physical,baseline=baseline,cut=cut)


def cap(m,depth,ones):
    zeros=depth-ones
    max_x=m-zeros
    max_x-=max_x%2
    return max_x+(m-max_x)//2


def bound(xs,entries,m,cut):
    j=len(xs)
    return sum(w+cap(m,j,sum((a&x).bit_count()%2 for x in xs))<cut for a,w in entries)


def payload(model,data,reps,light=None):
    text=f'{model.m} {len(data["labels"])} {model.d}\n'
    text+=' '.join(map(str,data['distances']))+'\n'
    text+=' '.join(map(str,data['Gamma']))+'\n'
    text+=str(len(reps))+'\n'+' '.join(map(str,reps))+'\n'
    if light is not None:
        text+=str(len(light['entries']))+'\n'
        text+='\n'.join(f'{a} {w}' for a,w in light['entries'])+'\n'
    return text
