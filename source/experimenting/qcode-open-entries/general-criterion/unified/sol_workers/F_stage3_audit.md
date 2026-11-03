**No correctness bugs found in the requested light-count prune or final equality test.** The audit was read-only and used short controls and evidence replays; conclusions remain conditional on the supplied suffix census and profile certificates.

1. **Budgets and the entire local stabilizer are correct.**

   | Frame | Profile | `(B6,B7,L)` |
   |---|---|---|
   | 688 | `(1,0)` | `(0,0,0)` |
   | 688 | `(1,1)` | `(0,1,1)` |
   | 688 | `(2,0)` | `(1,0,1)` |
   | 888 | `(0,0)` | `(0,0,0)` |
   | 888 | `(0,1)` | `(0,1,1)` |
   | 888 | `(0,2)` | `(0,2,2)` |
   | 888 | `(1,0)` | `(1,0,1)` |
   | 888 | `(1,1)` | `(1,1,2)` |
   | 888 | `(2,0)` | `(2,0,2)` |

   For 888, an extra local coset has four distinct words totaling weight 24. Minimum six forces four weight-six words, exceeding every budget. For 688, their total is at most 26, whereas at most one light word and minimum six require at least `6+3·8=30`. Thus the rank-two plane is the entire local stabilizer, including these minimum-six profiles. [PROOF.md:44](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/PROOF.md:44)

2. **Partial-depth counts are sound and do not double count.** `caps[j]` takes a maximum over every unknown prefix coordinate. Therefore `U<8` guarantees one light translate per distinct nonzero coefficient; `U<=6`, together with minimum six, guarantees weight six. The square nonsingular count activates each nonzero Z subset once and counts its four distinct translates. Its accumulated count and the cap count are tested separately, never summed. [enumerate_z.cpp:53](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/enumerate_z.cpp:53), [enumerate_z.cpp:104](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/enumerate_z.cpp:104)

3. **The full-Z combined count is a valid lower bound.** For `a∈Ann(X)`, expressing `Γa` in independent X determines z exactly, including `z=0` suffix-only words. Every translate is counted. `known[a]` excludes precisely those coefficients from the subsequent cap pass, so the contributions are disjoint. Excess exact weight-seven counts and combined guaranteed weight-six/light counts are legitimate exclusions. [enumerate_z.cpp:66](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/enumerate_z.cpp:66)

4. **Post-Y equality is correct.** `solve(...,find_all=True)` exhausts Y, retaining diagonal freedom. `graph_counts` uses the correct adjoint coordinates `(Yᵀa,Xᵀa)` and all four translates for every `a≠0`. Acceptance requires both the outside minimum and exact `(B6,B7)` equality. [common.py:55](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/common.py:55), [run.py:25](/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/general-criterion/unified/mixed_frames/stage3/run.py:25)

5. **Inherited 888 exclusions have the stated limited scope.**
   - `(0,1)` and `(0,2)`: mu7 excludes hull-zero classes **7,32,41,50,52**, plus all **41 singular classes**. Ten hull-zero classes remain open.
   - `(0,0)`: mu8 excludes all 15 hull-zero classes; the singular obstruction completes coverage of all 56.
   
   The singular argument is valid: the balanced kernel gives `hull=2q` and suffix-only stabilizer dimension q. Hull two forces a forbidden weight-six radical word. Hull six requires eight suffix-only words; their even weights and minimum seven force distance eight, contradicting the length-nine quaternary Plotkin bound `|N|≤6`. These exclusions do not transfer to minimum-six profiles.

Independent controls checked physical histograms, every depth of small known completions, singular nontrivial `Ann(X)` including suffix-only coefficients, insufficient-budget rejection, and final equality rejecting oversized budgets. Saved mu7/mu8 hashes, inputs, orbit coverage and terminal logs also replayed successfully.

**A rigorously safe symmetry reduction is available.**

- Use prefix **`S8 × C2`**: arbitrary site permutations and the uniform Clifford `X↔Y, Z→Z`. These preserve R, the pure-Z Lagrangian, every cost and every four-translate weight multiset. I checked all eight generating transformations against all 4,096 prefix classes.
- Use physical suffix automorphisms preserving C. If their coefficient action is T, their **syndrome action is `T⁻ᵀ`**, including singular Γ. I generated and checked these actions for all 56 classes.
- Canonicalize X through its affine eight-point set
  ` {0,X0,…,X5,X0⊕…⊕X5}`.
  Sort the seven nonzero points and require their largest point to be the dependent XOR. This removes `7!=5,040` labelings. Canonicalize over the eight anchor choices for another factor eight: **40,320** altogether. Complete every Y for the retained X; do not impose arbitrary GL(6,2) changes or discard diagonal freedoms.

The measured first-image orbit counts give:

| 888 scheduling | Partitions |
|---|---:|
| Current labeled manifest | 497,720 |
| Suffix first-image orbits | **171,580** |
| Also apply inherited exclusions and mu7 first-image eligibility | **129,056** |

Those remaining roots can be grouped into **188 case jobs**, with explicit orbit-coverage records. Prefix canonicalization separately removes a factor **40,320** from fully labeled injective X maps. Dividing 497,720 by 40,320 would **not** estimate the number of root partitions, and the available pilots do not support an exhaustive runtime forecast.