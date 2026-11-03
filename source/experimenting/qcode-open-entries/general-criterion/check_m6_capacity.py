"""Independent exploratory replay of the m=6 hull/capacity obstruction.
Uses BFS for coset distances, and independently builds all eligible prefix kernels.
Does not modify either submitted manuscript or the release package.
"""
from pathlib import Path
from collections import Counter, deque
from itertools import combinations, product
import json

def span(basis):
    out = [0]
    for x in basis:
        out += [y ^ x for y in out]
    return out

def nullspace(rows, width):
    pivots = {}
    for x in rows:
        for b in sorted(pivots, reverse=True):
            if (x >> b) & 1:
                x ^= pivots[b]
        if x:
            b = x.bit_length() - 1
            for c in list(pivots):
                if (pivots[c] >> b) & 1:
                    pivots[c] ^= x
            pivots[b] = x
    return [sum([1 << j] + [1 << b for b, row in pivots.items() if (row >> j) & 1])
            for j in range(width) if j not in pivots]

def weight(x, length):
    return ((x & ((1 << length)-1)) | (x >> length)).bit_count()

def symp(x, y, length):
    mask = (1 << length)-1
    return (((x & mask) & (y >> length)) ^ ((y & mask) & (x >> length))).bit_count() & 1

m = 6
p = 63 << 6
canon = lambda x: min(x, x ^ p)
E = sorted({canon(x) for x in range(4096) if (x & 63).bit_count() % 2 == 0})
cost = {x: min(weight(x,6),weight(x^p,6)) for x in E}
high = [x for x in E if cost[x] >= 5]
highset = set(high)
planes = {tuple(sorted((0,x,y,canon(x^y)))) for x,y in combinations(high,2) if canon(x^y) in highset}
# Every 3-dimensional space contains a plane. Check every possible extension.
for plane in planes:
    for x in high:
        if x not in plane:
            assert not all(canon(x^y) in highset for y in plane), 'A dimension-three kernel exists'
print(f'Prefix kernels: {len(high)} lines, {len(planes)} planes, no dimension-three subspace.')
profiles = {}
for K in [(0,)] + [(0,x) for x in high] + list(planes):
    u = (len(K)).bit_length()-1
    pairing_rank = 0 if all(symp(x,y,6)==0 for x in K for y in K) else 2
    seen, hist = set(), Counter()
    for e in E:
        if e in seen:
            continue
        block = {canon(e^x) for x in K}
        seen.update(block)
        minimum = min(cost[x] for x in block)
        if minimum:
            hist[minimum] += 1
    assert sum(hist.values()) == 1024//len(K)-1
    profiles.setdefault((u,pairing_rank),set()).add(tuple(sorted(hist.items())))
print('Distinct quotient-cost profiles:', {str(k):len(v) for k,v in profiles.items()})
root = Path(__file__).resolve().parents[3]
classes = json.loads((root/'quantum-code-bounds/verification/s7_37_classes.json').read_text())
results = []
for number, generators in enumerate(classes, 1):
    assert len(set(span(generators))) == 1024
    assert min(weight(x,10) for x in span(generators) if x) >= 5
    swapped = [(x>>10) | ((x&1023)<<10) for x in generators]
    P = nullspace(swapped,20)
    assert len(P) == 10
    assert all(symp(x,y,10)==0 for x in P for y in generators)
    gram = [sum(symp(x,y,10)<<j for j,y in enumerate(P)) for x in P]
    hull = span(nullspace(gram,10))
    h = len(hull).bit_length()-1
    def syndrome(x):
        return sum(symp(x,y,10)<<j for j,y in enumerate(P))
    steps = list({syndrome(a<<site) for site in range(10) for a in (1,1<<10,1|(1<<10))})
    distances = [-1]*1024
    distances[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for step in steps:
            y = x ^ step
            if distances[y] == -1:
                distances[y] = distances[x]+1
                queue.append(y)
    assert min(distances) >= 0
    survives = False
    for u in range(3):
        # Here dim(E)=dim(P)=10, so dim(K)=u and rank(K)=2u-h.
        candidate_profiles = profiles.get((u,2*u-h),set())
        if not candidate_profiles or u > h:
            continue
        if u == 0:
            kernels = [(0,)]
        elif u == 1:
            kernels = [(0,x) for x in hull if x]
        else:
            kernels = {tuple(sorted((0,x,y,x^y))) for x,y in combinations(hull,2) if x and y and x!=y}
        for N in kernels:
            hist = Counter(distances[x] for x in range(1,1024)
                           if all((x & y).bit_count()%2==0 for y in N))
            assert sum(hist.values()) == 1024//(1<<u)-1
            for profile in candidate_profiles:
                prefix_hist = dict(profile)
                if all(sum(count for value,count in hist.items() if value<=j)
                       <= sum(count for value,count in prefix_hist.items() if value>=5-j)
                       for j in range(5)):
                    survives = True
    assert not survives, f'Class {number} survives'
    results.append({'class':number,'hull':h,'verdict':'excluded by hull/capacity'})
print(f'PASS: all {len(results)} classes excluded by the necessary hull/capacity conditions.')

# Positive control: the known [[16,4,5]] stabilizer must pass the same inequality.
rows = [line.split('|') for line in (Path(__file__).resolve().parents[1]/'shorten/q1645.txt').read_text().splitlines() if line.strip()]
Sbasis = [sum(int(x)<<i for i,x in enumerate(a.split())) |
          (sum(int(x)<<i for i,x in enumerate(b.split()))<<16) for a,b in rows]
Swords = span(Sbasis)
assert len(set(Swords)) == 4096
assert all(symp(x,y,16)==0 for x in Sbasis for y in Sbasis)
pword = min((x for x in Swords if x),key=lambda x:weight(x,16))
assert weight(pword,16) == 6
T = [i for i in range(16) if (((pword&65535)|(pword>>16))>>i)&1]
W = [i for i in range(16) if i not in T]
def restrict(x, sites):
    return sum(((x>>i)&1)<<j for j,i in enumerate(sites)) | (sum(((x>>(16+i))&1)<<j for j,i in enumerate(sites))<<len(sites))
def dual(basis, length):
    mask = (1<<length)-1
    return nullspace([(x>>length)|((x&mask)<<length) for x in basis],2*length)
# Explicitly verify no error of weight below five commutes with S.
for size in range(1,5):
    for sites in combinations(range(16),size):
        for letters in product((1,1<<16,1|(1<<16)),repeat=size):
            error = sum(letter<<site for letter,site in zip(letters,sites))
            assert any(symp(error,s,16) for s in Sbasis)
Cbasis = dual([restrict(s,W) for s in Sbasis],10)
assert len(Cbasis)==9 and min(weight(x,10) for x in span(Cbasis) if x)>=5
Pbasis = dual(Cbasis,10)
assert len(Pbasis)==11
Pcoordinates = {0:0}
for j,x in enumerate(Pbasis):
    Pcoordinates.update({word^x:bits|(1<<j) for word,bits in list(Pcoordinates.items())})
N = {restrict(s,W) for s in Swords if restrict(s,T)==0}
Ncoeff = [Pcoordinates[x] for x in N]
u = len(N).bit_length()-1
Lbasis = dual([restrict(s,T) for s in Sbasis],6)
L = span(Lbasis)
assert restrict(pword,T) in L
Rword = restrict(pword,T)
allprefix = [x for x in range(4096) if symp(x,Rword,6)==0]
seen, prefix_hist = set(), Counter()
for x in allprefix:
    if x in seen: continue
    block = {x^y for y in L}
    seen.update(block)
    minimum = min(weight(y,6) for y in block)
    if minimum: prefix_hist[minimum]+=1
codim = len(Pbasis)
def control_syn(x):
    return sum(symp(x,y,10)<<j for j,y in enumerate(Pbasis))
steps = {control_syn(letter<<site) for site in range(10) for letter in (1,1<<10,1|(1<<10))}
distances = [-1]*(1<<codim); distances[0]=0; queue=deque([0])
while queue:
    x=queue.popleft()
    for step in steps:
        y=x^step
        if distances[y]<0:
            distances[y]=distances[x]+1;queue.append(y)
assert min(distances)>=0
suffix_hist=Counter(distances[x] for x in range(1,1<<codim) if all((x&y).bit_count()%2==0 for y in Ncoeff))
assert sum(prefix_hist.values())==sum(suffix_hist.values())
assert all(sum(count for value,count in suffix_hist.items() if value<=j)
           <=sum(count for value,count in prefix_hist.items() if value>=5-j) for j in range(5))
print(f'PASS: known [[16,4,5]] positive control satisfies the capacity inequalities (suffix kernel dimension {u}).')
