"""Direct labeled Python reference with optional prefix/depth and graph prunes.

The mu=0 traversal uses the stage-1 predicates verbatim. No bitsets or C++
linear algebra are used. Depth-limited leaves denote partial frontiers only.
"""
import time
from mixed_frames.enumerate_z import restricted_rank
from .common import graph_partial_ok, graph_z_ok, nonsingular_graph_weights


def enumerate_reference(model, frame, *, mu=0, prefix=(), depth=None,
                        node_limit=None, leaf_limit=None, use_rank=True, use_domains=True):
    p, Q, N = frame.p, 1 << frame.p, 1 << model.b
    depth = p if depth is None else depth
    sectors = [s for s in range(1,Q) if s.bit_count() <= 2]
    initial = {s: [a for a in range(N) if model.distances[a]+frame.costs[s << p] >= model.d]
               for s in sectors}
    leaves, nodes = [], 0
    reason = None
    graph_weights = nonsingular_graph_weights(model) if mu else None
    class Limit(Exception):
        pass
    def dfs(xs, images, domains):
        nonlocal nodes
        if node_limit is not None and nodes >= node_limit:
            raise Limit('nodes')
        nodes += 1
        j = len(xs)
        if j == depth:
            if j == p and not graph_z_ok(model,frame,xs,mu):
                return
            leaves.append(xs[:])
            if leaf_limit is not None and len(leaves) >= leaf_limit:
                raise Limit('leaves')
            return
        for beta in ([prefix[j]] if j < len(prefix) else range(N)):
            if beta in images:
                continue
            added = [a ^ beta for a in images]
            if any(frame.costs[(1<<j)|z]+model.distances[a]<model.d for z,a in enumerate(added)):
                continue
            if graph_weights is not None and any(frame.costs[(1<<j)|z]+graph_weights[a]<mu for z,a in enumerate(added)):
                continue
            new_xs = xs+[beta]
            if use_rank and restricted_rank(new_xs,model.Gamma)>2*(p-j-1):
                continue
            if not graph_partial_ok(model,frame,new_xs,mu):
                continue
            nd = domains
            if use_domains:
                nd = {}
                for s in sectors:
                    nd[s] = [candidate for candidate in domains[s]
                             if all(model.distances[candidate^a]+frame.costs[(s<<p)|(1<<j)|z]>=model.d
                                    for z,a in enumerate(added))]
                    if not nd[s]:
                        break
                if len(nd) != len(sectors) or any(not v for v in nd.values()):
                    continue
            dfs(new_xs, images+added, nd)
    try:
        if node_limit == 0:
            raise Limit('nodes')
        if leaf_limit == 0:
            raise Limit('leaves')
        if (not use_domains or all(initial.values())) and graph_partial_ok(model,frame,[],mu):
            dfs([], [0], initial if use_domains else {})
        status = 'COMPLETE'
    except Limit as exc:
        status, reason = 'LIMIT', str(exc)
    return dict(status=status, complete=status == 'COMPLETE', nodes=nodes,
                leaf_count=len(leaves), leaves=leaves, reason=reason)
