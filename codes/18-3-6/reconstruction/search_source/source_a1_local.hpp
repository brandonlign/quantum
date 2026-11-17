#pragma once
#include <bits/stdc++.h>
using namespace std;
inline int physp(int a,int b,int n){int m=(1<<n)-1;return __builtin_parity((a&m)&(b>>n))^__builtin_parity((b&m)&(a>>n));}
inline int wt(int a,int n){return __builtin_popcount((a&((1<<n)-1))|(a>>n));}
struct Basis{int b[40]{},labels[40]{};bool add(int x,int lab=0){while(x){int k=31-__builtin_clz(x);if(b[k]){x^=b[k];lab^=labels[k];}else{b[k]=x;labels[k]=lab;return true;}}return false;}int coords(int x)const{int z=0;while(x){int k=31-__builtin_clz(x);assert(b[k]);x^=b[k];z^=labels[k];}return z;}};
inline vector<int> span(const vector<int>&r){vector<int>w(1<<r.size());for(int i=1;i<(int)w.size();i++)w[i]=w[i&(i-1)]^r[__builtin_ctz(i)];return w;}
inline vector<int> ker(vector<int>eq,int width){Basis b;for(int x:eq)b.add(x);vector<int>out;for(int j=0;j<width;j++)if(!b.b[j]){int x=1<<j;for(int k=0;k<width;k++)if(b.b[k]&&__builtin_parity(b.b[k]&x))x^=1<<k;out.push_back(x);}return out;}
inline vector<int> perp(const vector<int>&r,int n){vector<int>eq;for(int x:r)eq.push_back((x>>n)|((x&((1<<n)-1))<<n));return ker(eq,2*n);}
inline vector<int> hyper(const vector<int>&c,int n){auto V=perp(c,n);vector<int>H;while(true){int ii=-1,jj=-1;for(int i=0;i<(int)V.size()&&ii<0;i++)for(int j=0;j<i;j++)if(physp(V[i],V[j],n)){ii=i;jj=j;break;}if(ii<0)return H;int a=V[ii],b=V[jj];vector<int>W;for(int k=0;k<(int)V.size();k++)if(k!=ii&&k!=jj){int v=V[k];W.push_back(v^(physp(v,b,n)?a:0)^(physp(v,a,n)?b:0));}H.push_back(a);H.push_back(b);V=move(W);}}
struct Local{int id,n,k,s;vector<int>C,Q;array<uint8_t,2048>cost;};
inline Local prepare(int id,int n,int s,const vector<int>&C){Local o;o.id=id;o.n=n;o.s=s;o.k=C.size();o.C=C;o.Q=hyper(C,n);assert(o.Q.size()==8);Basis b;for(int v:C)assert(b.add(v));for(int v:o.Q)assert(b.add(v));vector<int>A(C.begin(),C.begin()+s);for(int v:perp(A,n))if(b.add(v))o.Q.push_back(v);assert(o.Q.size()==11);auto cw=span(C),qw=span(o.Q);o.cost.fill(255);for(int v:cw)if(v)assert(wt(v,n)>=6);for(int i=0;i<2048;i++){int best=255;for(int c:cw)best=min(best,wt(qw[i]^c,n));o.cost[i]=best;assert((i==0)==(best==0));}return o;}
inline void put32(ostream&f,uint32_t x){f.write(reinterpret_cast<const char*>(&x),4);}
inline uint32_t get32(istream&f){uint32_t x;f.read(reinterpret_cast<char*>(&x),4);assert(f);return x;}
inline void write_local(ostream&f,const Local&o){put32(f,o.id);for(int x:o.C)put32(f,x);for(int x:o.Q)put32(f,x);f.write(reinterpret_cast<const char*>(o.cost.data()),2048);}
inline Local read_local(istream&f,int n,int k,int s){Local o;o.id=get32(f);o.n=n;o.k=k;o.s=s;o.C.resize(k);o.Q.resize(11);for(int&x:o.C)x=get32(f);for(int&x:o.Q)x=get32(f);f.read(reinterpret_cast<char*>(o.cost.data()),2048);assert(f);return o;}
inline void json_vec(ostream&o,const vector<int>&v){o<<"[";for(int i=0;i<(int)v.size();i++)o<<(i?",":"")<<v[i];o<<"]";}
inline void json_local(ostream&f,const Local&o){f<<"{\"id\":"<<o.id<<",\"n\":"<<o.n<<",\"s\":"<<o.s<<",\"C\":";json_vec(f,o.C);f<<",\"Q\":";json_vec(f,o.Q);f<<",\"cost\":";json_vec(f,vector<int>(o.cost.begin(),o.cost.end()));f<<"}";}
