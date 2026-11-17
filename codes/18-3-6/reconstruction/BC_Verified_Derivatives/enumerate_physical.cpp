#include <cstdint>
#include <fstream>
#include <iostream>
#include <set>
#include <stdexcept>
using namespace std;using U=uint64_t;
int wt(U x){return __builtin_popcountll((x&((1ULL<<18)-1))|(x>>18));}
int main(int argc,char**argv){try{ifstream f(argv[1]);U S[15],N[21];for(U&x:S)f>>x;for(U&x:N)f>>x;if(!f)throw runtime_error("input");set<U>stab;U x=0;stab.insert(0);for(unsigned j=1;j<(1u<<15);j++){x^=S[__builtin_ctz(j)];stab.insert(x);}if(stab.size()!=1u<<15)throw runtime_error("S rank");uint64_t A[19]{},B[19]{},L[19]{};U witness[19]{},normalw[19]{};for(U s:stab)A[wt(s)]++;x=0;B[0]++;for(unsigned j=1;j<(1u<<21);j++){x^=N[__builtin_ctz(j)];int w=wt(x);B[w]++;normalw[w]=x;if(!stab.count(x)){L[w]++;witness[w]=x;}}cout<<"weight stabilizer normalizer logical logical_witness\n";int dn=99,dl=99;for(int w=0;w<=18;w++){cout<<w<<' '<<A[w]<<' '<<B[w]<<' '<<L[w]<<' '<<witness[w]<<'\n';if(w&&B[w])dn=min(dn,w);if(L[w])dl=min(dl,w);if(B[w]!=A[w]+L[w])throw runtime_error("enumerators");}cout<<"PASS all 2097152 normalizer words; minimum "<<dn<<" logical distance "<<dl<<" witness "<<witness[dl]<<endl;return 0;}catch(exception&e){cerr<<e.what()<<endl;return 1;}}
