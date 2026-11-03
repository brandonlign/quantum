"""Physical prefix actions and finite matrix closure with full lifts."""
from model import coordinates,lift_coordinate,span,weight,gram,reduced_basis
from m8r1.suffix_symmetries import automorphisms


def closure(generators,width):
    identity=tuple(1<<i for i in range(width));seen={identity};todo=[identity]
    while todo:
        cur=todo.pop()
        for g in generators:
            new=tuple(lift_coordinate(a,cur) for a in g)
            if new not in seen:seen.add(new);todo.append(new)
    return sorted(seen)


def prefix_group(frame):
    p=frame.p;r=len(frame.R);fullword=((1<<frame.m)-1)<<frame.m
    maps,order=automorphisms(frame.R,frame.m,[fullword]);generators=[]
    for apply in maps:
        U=[coordinates(apply(w),frame.R+frame.E)>>r for w in frame.E]
        assert all(v<1<<p for v in U[:p])
        assert len(reduced_basis(U))==2*p
        for e in range(1<<(2*p)):
            w=lift_coordinate(e,frame.E);v=lift_coordinate(lift_coordinate(e,U),frame.E)
            assert sorted(weight(w^h,frame.m) for h in span(frame.R))==sorted(weight(v^h,frame.m) for h in span(frame.R))
        generators.append(U)
    identity=tuple(1<<i for i in range(2*p));lifts={identity[:p]:identity};todo=[identity]
    while todo:
        cur=todo.pop()
        for g in generators:
            new=tuple(lift_coordinate(a,cur) for a in g)
            if new[:p] not in lifts:lifts[new[:p]]=new;todo.append(new)
    return dict(physical_group_order=order,full_generators=generators,
                Z_group=[list(g) for g in sorted(lifts)],full_lifts=[list(lifts[g]) for g in sorted(lifts)])
