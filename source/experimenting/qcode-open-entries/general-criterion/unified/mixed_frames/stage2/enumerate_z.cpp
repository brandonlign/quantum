#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
using Bits = array<uint64_t,64>;
int p,b,Q,N,d,mu,rank_on,domains_on,depth;
bool ack=false;
vector<int> cost,dist,gamma,wt,xs,prefix,sectors;
vector<int> graph_weight;
vector<vector<int>> caps,light;
vector<Bits> good;
uint64_t nodes=0,leaves=0,tests=0,rank_prunes=0,sector_prunes=0,mu_prunes=0;
uint64_t node_limit=UINT64_MAX,leaf_limit=UINT64_MAX;
double seconds=-1;
string reason;
auto started=chrono::steady_clock::now();
int parity(int v){return __builtin_parity((unsigned)v);}
double elapsed(){return chrono::duration<double>(chrono::steady_clock::now()-started).count();}
void limits(){if(seconds>=0 && elapsed()>=seconds) throw string("time");}
vector<int> annihilate(const vector<int>& ann,int beta){
    vector<int> out; int pivot=0;
    for(int a:ann) if(parity(a&beta)){pivot=a;break;}
    for(int a:ann) if(a!=pivot) out.push_back(a ^ (parity(a&beta)?pivot:0));
    return out;
}
int lift(int a,const vector<int>& v){int out=0; for(unsigned i=0;i<v.size();i++) if(a>>i&1) out^=v[i]; return out;}
int matrix_rank(vector<int> rows){int piv[24]={},r=0;for(int a:rows){while(a){int j=31-__builtin_clz((unsigned)a);if(piv[j])a^=piv[j];else{piv[j]=a;r++;break;}}}return r;}
int restricted(const vector<int>& ann){vector<int> rows; for(int a:ann){int row=0;for(unsigned j=0;j<ann.size();j++)row|=parity(a&lift(ann[j],gamma))<<j;rows.push_back(row);}return matrix_rank(rows);}
uint64_t permute(uint64_t v,int s){
    if(s&1) v=((v&0x5555555555555555ULL)<<1)|((v>>1)&0x5555555555555555ULL);
    if(s&2) v=((v&0x3333333333333333ULL)<<2)|((v>>2)&0x3333333333333333ULL);
    if(s&4) v=((v&0x0f0f0f0f0f0f0f0fULL)<<4)|((v>>4)&0x0f0f0f0f0f0f0f0fULL);
    if(s&8) v=((v&0x00ff00ff00ff00ffULL)<<8)|((v>>8)&0x00ff00ff00ff00ffULL);
    if(s&16)v=((v&0x0000ffff0000ffffULL)<<16)|((v>>16)&0x0000ffff0000ffffULL);
    if(s&32)v=(v<<32)|(v>>32);
    return v;
}
bool intersect(Bits& dst,int need,int shift){
    if(need<=0) return true;
    bool any=false; int blocks=(N+63)/64;
    for(int k=0;k<blocks;k++){dst[k]&=permute(good[need][k^(shift>>6)],shift&63);any|=dst[k]!=0;}
    return any;
}
bool graph_partial(int j){
    if(mu<=0) return true;
    for(int a:light[j]){int x=0;for(int i=0;i<j;i++)x|=parity(xs[i]&a)<<i;
        if(caps[j][x]+wt[a]<mu) return false;}
    return true;
}
bool graph_z(const vector<int>& ann){
    if(mu<=0) return true;
    int piv[24]={},coord[24]={};
    for(int i=0;i<p;i++){int v=xs[i],c=1<<i;while(v){int j=31-__builtin_clz((unsigned)v);if(piv[j]){v^=piv[j];c^=coord[j];}else{piv[j]=v;coord[j]=c;break;}}}
    vector<int> as={0};for(int g:ann){int n=as.size();for(int i=0;i<n;i++)as.push_back(as[i]^g);}
    for(int a:as) if(a){int v=lift(a,gamma),z=0;while(v){int j=31-__builtin_clz((unsigned)v);if(!piv[j]) return false;v^=piv[j];z^=coord[j];}
        // a!=0 labels a nonzero physical suffix word, including z==0.
        if(cost[z]+wt[a]<mu) return false;}
    return true;
}
void dfs(const vector<int>& images,const vector<int>& ann,const vector<Bits>& domains){
    limits();if(nodes>=node_limit)throw string("nodes");nodes++;
    int j=xs.size();
    if(j==depth){
        if(j==p && !graph_z(ann)){mu_prunes++;return;}
        leaves++;cout<<"{\"type\":\"leaf\",\"X\":[";for(int i=0;i<j;i++)cout<<(i?",":"")<<xs[i];cout<<"]}"<<endl;
        if(ack){string line;if(!getline(cin,line)||line!="OK")throw string("completion_pipeline");}
        if(leaves>=leaf_limit)throw string("leaves");return;
    }
    int lo=0,hi=N;if(j<(int)prefix.size()){lo=prefix[j];hi=lo+1;}
    // All newly activated pure-Z tests are translated bitset intersections.
    Bits options=good[0];
    for(unsigned z=0;z<images.size();z++)if(!intersect(options,d-cost[(1<<j)|z],images[z]))return;
    for(int a:images)options[a/64]&=~(1ULL<<(a%64));
    for(int block=lo/64;block<=(hi-1)/64;block++){
      uint64_t choices=options[block];
      while(choices){int bit=__builtin_ctzll(choices);choices&=choices-1;int beta=64*block+bit;if(beta<lo||beta>=hi)continue;
        tests++;if((tests&4095)==0)limits();
        vector<int> added;bool ok=true;
        for(unsigned z=0;z<images.size();z++){int a=images[z]^beta;added.push_back(a);
            if(!graph_weight.empty() && cost[(1<<j)|z]+graph_weight[a]<mu){ok=false;break;}}
        if(!ok){mu_prunes++;continue;}
        if(!ok)continue;
        auto newann=annihilate(ann,beta);
        if(rank_on && restricted(newann)>2*(p-j-1)){rank_prunes++;continue;}
        xs.push_back(beta);
        if(!graph_partial(j+1)){mu_prunes++;xs.pop_back();continue;}
        auto nd=domains;
        if(domains_on)for(unsigned k=0;k<sectors.size() && ok;k++){
            int s=sectors[k];
            for(unsigned z=0;z<added.size();z++)if(!intersect(nd[k],d-cost[(s<<p)|(1<<j)|z],added[z])){ok=false;break;}
        }
        if(!ok){sector_prunes++;xs.pop_back();continue;}
        auto ni=images;ni.insert(ni.end(),added.begin(),added.end());
        dfs(ni,newann,nd);xs.pop_back();
      }
    }
}
int main(int argc,char**argv){
    if(argc<2){cerr<<"input required\n";return 2;}
    ifstream in(argv[1]);in>>p>>b>>d>>mu>>rank_on>>domains_on;
    if(!in||p<1||p>6||b<1||b>12){cerr<<"invalid input\n";return 2;}
    Q=1<<p;N=1<<b;depth=p;
    cost.resize(Q*Q);dist.resize(N);gamma.resize(b);wt.resize(N);
    for(int&v:cost)in>>v;for(int&v:dist)in>>v;for(int&v:gamma)in>>v;for(int&v:wt)in>>v;
    if(!in){cerr<<"truncated input\n";return 2;}
    for(int k=2;k<argc;k+=2){if(k+1==argc)return 2;string a=argv[k],v=argv[k+1];
        if(a=="--nodes")node_limit=stoull(v);else if(a=="--leaves")leaf_limit=stoull(v);else if(a=="--seconds")seconds=stod(v);
        else if(a=="--depth")depth=stoi(v);else if(a=="--prefix"){stringstream ss(v);string t;while(getline(ss,t,','))prefix.push_back(stoi(t));}
        else if(a=="--ack")ack=stoi(v)!=0;
        else return 2;
    }
    if(depth<0||depth>p||(int)prefix.size()>depth||any_of(prefix.begin(),prefix.end(),[](int a){return a<0||a>=N;}))return 2;
    good.resize(32);for(int need=0;need<32;need++)for(int a=0;a<N;a++)if(dist[a]>=need)good[need][a/64]|=1ULL<<(a%64);
    // Only for nonsingular Gamma: compatibility makes M invertible, and the
    // unique stabilizer suffix of pure-Z e is P Gamma^{-1} M e. This stronger
    // early test is never applied to singular syndrome forms.
    if(mu>0 && b==2*p && matrix_rank(gamma)==b){
        graph_weight.assign(N,0);for(int a=0;a<N;a++)graph_weight[lift(a,gamma)]=wt[a];
    }
    caps.resize(p+1);light.resize(p+1);
    for(int j=0;j<=p;j++){caps[j].assign(1<<j,0);for(int x=0;x<Q;x++)for(int z=0;z<Q;z++)caps[j][x&((1<<j)-1)]=max(caps[j][x&((1<<j)-1)],cost[(x<<p)|z]);
        int mincap=*min_element(caps[j].begin(),caps[j].end());for(int a=1;a<N;a++)if(wt[a]+mincap<mu)light[j].push_back(a);}
    vector<Bits> domains;bool empty=false;
    for(int s=1;s<Q;s++)if(__builtin_popcount((unsigned)s)<=2){sectors.push_back(s);domains.push_back(good[max(0,d-cost[s<<p])]);if(all_of(domains.back().begin(),domains.back().end(),[](uint64_t a){return !a;}))empty=true;}
    vector<int> ann;for(int i=0;i<b;i++)ann.push_back(1<<i);
    string status="COMPLETE";
    try{limits();if(node_limit==0)throw string("nodes");if(leaf_limit==0)throw string("leaves");if(!domains_on||!empty){if(graph_partial(0))dfs({0},ann,domains);else mu_prunes++;}}
    catch(string r){status="LIMIT";reason=r;}
    cout<<"{\"type\":\"summary\",\"status\":\""<<status<<"\",\"complete\":"<<(status=="COMPLETE"?"true":"false")<<",\"reason\":\""<<reason<<"\",\"nodes\":"<<nodes<<",\"leaf_count\":"<<leaves<<",\"candidates\":"<<tests<<",\"rank_prunes\":"<<rank_prunes<<",\"sector_prunes\":"<<sector_prunes<<",\"mu_prunes\":"<<mu_prunes<<",\"seconds\":"<<elapsed()<<"}"<<endl;
}
