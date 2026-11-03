"""Compose compact singleton columns with optional exact MITM feasibility."""
import ast
from pathlib import Path
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
tree=ast.parse((ROOT/'fixed_kernel_mitm4/patch.py').read_text())
helper=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='helper' for t in n.targets))
def hrep(a,b):
 global helper
 assert helper.count(a)==1,a[:100];helper=helper.replace(a,b)
hrep('static uint64_t mdom[4][MAX_WORDS], mbase[MAX_WORDS], mH2[MAX_WORDS], minW[MAX_WORDS];',
'''static uint64_t mdom[4][MAX_WORDS], mbase[MAX_WORDS];
/* Pair sums need only delta>=d-2 and MUST NOT use compact singleton columns. */
static unsigned member_epoch[MAX_N],pair_epoch[MAX_N],node_epoch;
static unsigned char pair_good[MAX_N]; static int mspan;
static int inW(int v){return member_epoch[v]==node_epoch;}
static int pair_allowed(int v){
 if(pair_epoch[v]==node_epoch)return pair_good[v];
 int ok=1;
 for(int j=0;j<mspan;j++){
  int k=sizes[4][j],t=k+2<m-k-2?k+2:m-k-2;
  if(distance-t>0 && distances[v^sums[4][j]]<distance-t){ok=0;break;}
 }
 pair_epoch[v]=node_epoch;pair_good[v]=(unsigned char)ok;return ok;
}''')
helper=helper.replace('mbit(minW,','inW(').replace('mbit(mH2,','pair_allowed(')
hrep('if(!mbit(mdom[3],d))continue;',
     'int dc=label_column[d];if(dc<0 || !mbit(mdom[3],dc))continue;')
hrep('''    for(int w=0;w<word_count;w++){mbase[w]=UINT64_MAX;mH2[w]=UINT64_MAX;minW[w]=0;}
    if(N<64){mbase[0]=mH2[0]=(UINT64_C(1)<<N)-1;}''',
'''    for(int w=0;w<word_count;w++)mbase[w]=UINT64_MAX;
    if(ZN&63)mbase[word_count-1]=(UINT64_C(1)<<(ZN&63))-1;
    if(++node_epoch==0){memset(member_epoch,0,sizeof member_epoch);memset(pair_epoch,0,sizeof pair_epoch);node_epoch=1;}
    mspan=span_size;''')
hrep('''        if(distance-t2>0){const uint64_t *row=translated_row(t2,sums[depth][j]);for(int w=0;w<word_count;w++)mH2[w]&=row[w];}
        minW[sums[depth][j]>>6]|=UINT64_C(1)<<(sums[depth][j]&63);''',
'''        (void)t2;
        int v=sums[depth][j];member_epoch[v]=node_epoch;
        int col=label_column[v];if(col>=0)mbase[col>>6]&=~(UINT64_C(1)<<(col&63));''')
hrep('    for(int w=0;w<word_count;w++)mbase[w]&=~minW[w];','')
hrep('''        if(f<fixed_prefix_len){int label=fixed_prefix[f];uint64_t kept=out[label/64]&(UINT64_C(1)<<(label%64));
            for(int w=0;w<word_count;w++)out[w]=0;out[label/64]=kept;}''',
'''        if(f<fixed_prefix_len){int col=label_column[fixed_prefix[f]];if(col<0)return 0;
            uint64_t kept=out[col/64]&(UINT64_C(1)<<(col%64));
            for(int w=0;w<word_count;w++)out[w]=0;out[col/64]=kept;}''')
hrep('int v=64*w+__builtin_ctzll(b);','int v=column_label[64*w+__builtin_ctzll(b)];')
s=(ROOT/'kernel_z_compact/z_enum_compact.c').read_text()
for node in tree.body:
 if not isinstance(node,ast.Expr) or not isinstance(node.value,ast.Call) or not isinstance(node.value.func,ast.Name) or node.value.func.id!='rep':continue
 if isinstance(node.value.args[0],ast.Constant) and node.value.args[0].value.startswith('/* Exact streaming'):continue
 a,b=[eval(compile(ast.Expression(arg),'local-patch','eval'),{'helper':helper}) for arg in node.value.args]
 assert s.count(a)==1,a[:100];s=s.replace(a,b)
assert 'mH2' not in s and 'minW' not in s
(HERE/'z_enum_joint.c').write_text(s)
