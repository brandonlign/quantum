# Independent verification of the nine-site dimension-five distance-six cell

This verification uses an independently encoded physical equivalence graph and complete mass/lift checks. It replaces a second direct distance-six classification, provided the previously certified release16 seven-site dimension-four/five parent cells are accepted. It does not exclude a quantum code or establish novelty.

## Physical equivalence graph

For each additive code, create a vertex for every physical codeword, three letter vertices per site, and a separate site vertex joined to those three letters. Join a codeword vertex to its nonzero letter at each site. Give the three vertex kinds separate colors. A color-preserving automorphism is precisely a site permutation and a permutation of each site's three nonzero Pauli letters taking the code to itself. Each S3 letter permutation is binary-linear and symplectic on that site. Conversely every such physical equivalence acts on this graph. Codeword vertices have no extra kernel because their complete letter incidences are distinct. Thus its certificate classifies physical codes and its group order is the physical automorphism order. The earlier census instead used triangles without site vertices; the new representation is independent.

## Parent completeness

The certified release16 seven-site dimension-four/five cells have 43 and 1 classes. Their new graph certificate sets are compared with those used by the new lengthening check. The dimension-three seven-site cell is independently regenerated from every multiset of seven coordinate row spaces. There are 15 possible row spaces of dimension at most two in F2^3: zero, seven lines, seven planes. Their output maps are equivalent under local GL(2,2) exactly when the row spaces agree. Site permutations turn a code into a multiset of seven types. Every one of the 116280 multisets is tested; 2171 have dimension three and distance at least five, giving exactly the same 32 new graph certificates as the saved parent list.

## Eight-site labelled mass

Shortening any length-eight dimension-five distance-five code at the last site gives a length-seven parent of dimension 5-r, r=0,1,2. The extension is specified by an r-dimensional subspace of the quotient Pauli space by that parent, with every nonzero coset leader at least four. Conversely these data give a valid distance-five child.

The verifier computes all quotient minima by exhaustively assigning all 4^7 physical words to their symplectic syndrome labels. It checks the zero-syndrome kernel equals the entire parent code. For r=0 there is one extension. For r=1 there are three choices of nonzero last-site letter per deep coset. For r=2 there are six choices of the isomorphism from the last-site Pauli space to each deep two-space. Consequently the labelled r=2 extension count is the number of ordered deep pairs a,b whose sum is also deep. It is independently computed both by direct ordered-pair counting and by the Fourier identity

    count = (1/|Q|) sum_u Fourier(1_deep)(u)^3.

For each parent P, multiply this count by 7!6^7/|Aut(P)|. Summing gives the total number of labelled children, with the last site fixed. Each labelled child has a unique shortening parent and extension, so this is a complete universe count.

Every saved child is independently checked for rank five, distance five, new graph certificate distinctness, and automorphism order. The sum 8!6^8/|Aut(child)| equals the same full universe count:

    1593099450132480.

Because the child list consists of disjoint valid physical orbits whose total mass equals the full universe, it contains every child orbit.

## Nine-site puncture/lift coverage

Puncturing a length-nine dimension-five distance-six code at one site preserves dimension: a nonzero puncture kernel would contain a weight-one codeword. The punctured code has distance at least five and is physically equivalent to one of the now-checked eight-site parents. The new site's linear map is a pair of binary forms a,b on F2^5. All 1024 ordered pairs are tested. A pair is valid exactly when every weight-five coefficient word of the parent has a nonzero symbol at the new site.

Local GL(2,2) at the new site identifies pairs with the same row space: one zero pair, 31 lines, 155 planes, or 187 representatives. The verifier regenerates every valid canonical lift and checks its distance, hull and full nine-site graph certificate. Its certificate set and each individual saved class's multiplicity must match the source list, not merely the grand total or hull counts. Every source class must occur.

A PASS therefore covers all 299175 canonical lifts and all 1791031 unrestricted valid ordered pairs, and establishes the full 32865-class nine-site cell relative only to the certified release16 parent completeness. The independent direct distance-six census is no longer required as an additional coverage argument. An outside audit of this replacement remains appropriate before final publication.
