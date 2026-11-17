#include <array>
#include <cstdint>
#include <iostream>
#include <stdexcept>
using namespace std;using U=uint64_t;
const U S[15]={26277347457,59444428930,61476569286,11743039690,1976041682,19332759778,16496890048,63874044160,21506327040,66477130752,15987673088,20945838080,51785547776,16911630336,66846720};
unsigned syn[18][3];uint64_t counts[6]{},zeros[6]{};
int pairing(U a,U b){int p=0;for(int j=0;j<18;j++)p^=((a>>j)&1)*((b>>(18+j))&1)^((b>>j)&1)*((a>>(18+j))&1);return p;}
void go(int start,int remaining,int target,unsigned syndrome){if(!remaining){counts[target]++;if(!syndrome)zeros[target]++;return;}for(int j=start;j<=18-remaining;j++)for(int s=0;s<3;s++)go(j+1,remaining-1,target,syndrome^syn[j][s]);}
int main(){for(int i=0;i<18;i++)for(int t=1;t<=3;t++){U w=U(t&1)<<i|U(t>>1)<<(18+i);for(int k=0;k<15;k++)syn[i][t-1]|=pairing(w,S[k])<<k;}for(int w=1;w<=5;w++){go(0,w,w,0);cout<<"weight "<<w<<" errors "<<counts[w]<<" zero_syndrome "<<zeros[w]<<'\n';if(zeros[w])throw runtime_error("low-weight error");}U logical=34359939617ULL;for(U s:S)if(pairing(s,logical))throw runtime_error("logical witness syndrome");cout<<"PASS all physical errors of weights 1 through 5 have nonzero syndrome; witness has zero syndrome\n";}
