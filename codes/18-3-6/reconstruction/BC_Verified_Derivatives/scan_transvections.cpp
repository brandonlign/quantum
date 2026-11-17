#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
using namespace std;using U=uint64_t;constexpr U mask18=(1ULL<<18)-1;
int wt(U v){return __builtin_popcountll((v&mask18)|(v>>18));}int sp(U a,U b){return __builtin_parityll(((a&mask18)&(b>>18))^((b&mask18)&(a>>18)));}
vector<U>smallspan(vector<U>b,int limit){vector<U>out;U v=0;for(unsigned j=1;j<(1u<<b.size());j++){v^=b[__builtin_ctz(j)];if(wt(v)<=limit)out.push_back(v);}return out;}
int main(int argc,char**argv){ifstream f(argv[1]);vector<U>S(15),N(21);for(U&x:S)f>>x;for(U&x:N)f>>x;if(!f)return 1;bool pauli=argc>3;auto smallS=smallspan(S,11),smallN=smallspan(N,10);ofstream out(argv[2]);out<<"mask A7 D6 B_style C_style\n";long tested=0,eligible=0,valid=0,nb=0,nc=0;map<int,int>hist;
for(int a=0;a<18;a++)for(int b=a+1;b<18;b++)for(int c=b+1;c<18;c++)for(int d=c+1;d<18;d++){
 int sites[4]={a,b,c,d};for(int label=0;label<(pauli?81:1);label++){U v=0;int q=label;for(int i:sites){int p=pauli?1+q%3:2;q/=3;v|=U(p&1)<<i;v|=U(p>>1)<<(18+i);}tested++;bool ok=true;int a7=0;for(U x:smallS){int w=wt(x^(sp(x,v)?v:0));if(w<7){ok=false;break;}a7+=w==7;}if(!ok)continue;eligible++;int d6=0;for(U x:smallN){int w=wt(x^(sp(x,v)?v:0));if(w<6){ok=false;break;}d6+=w==6;}if(!ok)continue;valid++;hist[a7]++;bool B=a7==16&&d6==564,C=a7==12;nb+=B;nc+=C;out<<v<<' '<<a7<<' '<<d6<<' '<<B<<' '<<C<<'\n';if(B||C){cout<<"MATCH mask "<<v<<" A7 "<<a7<<" D6 "<<d6<<" B "<<B<<" C "<<C<<endl;}if(pauli&&nb&&nc){cout<<"BOUNDED_STOP after "<<tested<<" of247860 masks\n";goto done;}}}
done:out.flush();cout<<"DONE tested "<<tested<<" minS7 "<<eligible<<" valid_distance6 "<<valid<<" B_style "<<nb<<" C_style "<<nc<<" smallS "<<smallS.size()<<" smallN "<<smallN.size()<<'\n';for(auto [w,n]:hist)cout<<"A7 "<<w<<" valid_count "<<n<<'\n';}
