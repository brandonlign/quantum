from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'benchmarks/results/actual_17_screens.json').read_text())
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(10,4.5),layout='constrained')
for ax,b in zip(axes,data['branches']):
    values=[b['search_tasks_before_screens_after_symmetry'],b['search_tasks_after_screens']]
    bars=ax.bar(['Before tests','After tests'],values,color=['#9aa9b7','#225b85'],width=.57)
    ax.set_ylim(0,max(values)*1.23)
    ax.set_ylabel('Fixed-kernel search tasks')
    ax.set_title(('Hull-five branch' if b['branch'].endswith('h5') else 'Hull-three branch')+'\n'+f"{b['percent_search_tasks_avoided']:.1f}% fewer tasks")
    for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value+max(values)*.025,f'{value:,}',ha='center')
    ax.ticklabel_format(axis='y',style='plain')
fig.suptitle('Measured reductions in two seventeen-qubit proof branches',fontsize=14)
fig.savefig(ROOT/'benchmarks/results/seventeen_qubit_comparison.png',dpi=180)
fig.savefig(ROOT/'benchmarks/results/seventeen_qubit_comparison.svg')
print('Wrote comparison figures.')
