"""Hand-checkable congruence proving A4 even for [[16,5,>=5]] with A1=A2=A3=0.
Identities (all exact for any [[16,5,>=5]] code, A_w = #stabilizer words of weight w,
K_j = quaternary Krawtchouk polynomial of length 16, Sh_j = shadow counts, integers >= 0):
   Delta_0 = sum_w A_w                    = 2048
   Delta_j = sum_w K_j(w) A_w - 2048 A_j  = 0      (j = 1..4; S^perp has no words of wt 1..4 outside S)
   Gamma_j = sum_w (-1)^w K_j(w) A_w - 2048 Sh_j = 0
We find integers lam with L = sum lam_i * (identity_i) such that, modulo 2^e, every coefficient of
A_w (w != 4) and of Sh_j vanishes and the coefficient of A_4 is 2^(e-1).  Then
L's value  ≡ 2^(e-1) A_4 (mod 2^e) forces the parity of A_4.  Everything is re-verified exactly."""
from math import comb
from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp
n = 16
K = lambda j, w: sum((-1) ** a * 3 ** (j - a) * comb(w, a) * comb(n - w, j - a) for a in range(max(0, j - (n - w)), min(j, w) + 1))
# identity rows as dicts over symbols ('A',w) and ('Sh',j); value on the right
rows = [('Delta0', {('A', w): 1 for w in range(n + 1)}, 2048)]
for j in range(1, 5):
    r = {('A', w): K(j, w) for w in range(n + 1)}; r[('A', j)] -= 2048; rows.append((f'Delta{j}', r, 0))
for j in (0, 1, 3):
    r = {('A', w): (-1) ** w * K(j, w) for w in range(n + 1)}; r[('Sh', j)] = -2048; rows.append((f'Gamma{j}', r, 0))
# known: A0=1, A1=A2=A3=0 -> move to values; unknowns: A4..A16, Sh0, Sh1, Sh3
unk = [('A', w) for w in range(4, n + 1)] + [('Sh', j) for j in (0, 1, 3)]
B = Matrix([[r.get(u, 0) for u in unk] for _, r, _ in rows])
val = [v - r[('A', 0)] for _, r, v in rows]                    # value minus A0 contribution
best = None
for e in range(1, 30):
    Mod = 2 ** e
    # need lam: lam.B[:,col] ≡ 0 (cols != A4), lam.B[:,A4] ≡ Mod/2, mod Mod.  brute via SNF over Z of augmented
    target = [0] * len(unk); target[0] = Mod // 2
    # solve lam B ≡ target mod Mod  <=>  exists integer t: lam B - Mod t = target
    Aug = B.col_join(Mod * Matrix.eye(len(unk)))               # rows: lam (7) then t (len unk)
    D, U, V = smith_normal_decomp(Aug.T)                       # Aug.T * x = target
    c = U * Matrix(target); ok = True; sol = [0] * Aug.rows
    y = [0] * Aug.rows
    for i in range(D.rows):
        dii = D[i, i] if i < D.cols else 0
        if dii == 0:
            if c[i] != 0: ok = False; break
        else:
            if c[i] % dii: ok = False; break
            y[i] = c[i] // dii
    if ok:
        x = V * Matrix(y)
        lam = [int(v) % Mod for v in x[:len(rows)]]
        lam = [v - Mod if v > Mod // 2 else v for v in lam]
        best = (e, lam); break
e, lam = best; Mod = 2 ** e
coef = {u: sum(l * r.get(u, 0) for l, (_, r, _) in zip(lam, rows)) for u in [('A', w) for w in range(n + 1)] + [('Sh', j) for j in (0, 1, 3)]}
Lval = sum(l * v for l, (_, _, v) in zip(lam, rows))
print(f'modulus 2^{e} = {Mod}')
print('L =', ' '.join(f'{l:+d}*{nm}' for l, (nm, _, _) in zip(lam, rows) if l))
print('value of L:', Lval)
for u, c in coef.items(): print(f'  coef of {u[0]}_{u[1]}: {c}   mod {Mod}: {c % Mod}')
# exact check of the conclusion: with A0=1, A1=A2=A3=0:  L = coef_A0 + coef_A4*A4 + (multiples of Mod)
assert all(coef[u] % Mod == 0 for u in coef if u not in (('A', 0), ('A', 1), ('A', 2), ('A', 3), ('A', 4)))
assert coef[('A', 4)] % Mod == Mod // 2
rhs = (Lval - coef[('A', 0)]) % Mod
print(f'=> {Mod//2}*A4 ≡ {rhs} (mod {Mod})  =>  A4 ≡ {rhs // (Mod//2)} (mod 2)')
assert rhs in (0, Mod // 2)
