# Code A: a pure binary [[18,3,6]] stabilizer code

Code A has 18 physical qubits, encodes 3 logical qubits, and has distance exactly 6. Its 15 independent stabilizer generators commute. Every nonidentity element of its normalizer has weight at least 6, so the code is pure.

## Generators

Each line is an 18-qubit Pauli tensor product. Coordinates are written from qubit 1 to qubit 18. The generators may all be assigned positive Hermitian phases.

```text
XIIIZIIYZZZIIIIYZI
IYIZIIZYZIZIZZZIZY
IIXIZIZXIIZIYIIXZY
IZZYZZIYIZZZIZIIII
IIIIYZIXZIZZYIIYII
IZZIZXZXIIIIXZIZZX
IZIIZIYYZIZIYZZYIX
IIZZIIZZYZZIYZIYZZ
IIIZZZZIIXIIIIZXZI
IZZIZIIZIZYZZIZYZZ
IIZZZZIIIZZYIZZYII
IZZZZIIIIIIZZYIIZI
IZIZIZIZZZIIIIXXZZ
IIIIIIIIIIZZZZZZXX
ZZZZZZZZIIIIIIIIII
```

## Exact verification

- Binary generator rank: 15; all 105 distinct generator pairs commute.
- Normalizer rank: 21. All 2,097,152 normalizer vectors were enumerated exactly.
- Normalizer weights 1 through 5: zero.
- Weight-six normalizer vectors: 588, comprising 4 stabilizers and 584 nontrivial logical operators.
- An additional direct check finds nonzero syndrome for every one of the 2,353,347 physical Pauli errors of weights 1 through 5.
- All 19 coefficients of the quaternary MacWilliams transform of the stabilizer weight enumerator agree with the independently enumerated normalizer.

For example, XIIIIXIIIXIIXIIIXY is a weight-six logical operator: it commutes with every generator and is not in the stabilizer. This proves distance exactly 6, rather than merely a lower bound.

The stabilizer weight enumerator is

A(y) = 1 + 4y^6 + 130y^8 + 480y^9 + 1436y^10 + 2688y^11 + 4312y^12 + 6464y^13 + 7260y^14 + 5760y^15 + 2981y^16 + 992y^17 + 260y^18.

The normalizer weight enumerator is

B(y) = 1 + 588y^6 + 2144y^7 + 8538y^8 + 28960y^9 + 79436y^10 + 171968y^11 + 300808y^12 + 417344y^13 + 445700y^14 + 357856y^15 + 201501y^16 + 70304y^17 + 12004y^18.

The self-contained Python verifier `verify_standalone.py` uses only the standard library and no external files. Run it with Python 3, without the `-O` option. It reconstructs the normalizer from the displayed code, checks the full enumerators and MacWilliams identity, and verifies the explicit logical witness. The separate C++ program `check_errors.cpp` directly checks all errors of weights at most five.
