from pathlib import Path
import csv,json

ROOT=Path(__file__).resolve().parents[1]
def main():
    small=json.loads((ROOT/'benchmarks/results/framework_comparison.json').read_text())
    large=json.loads((ROOT/'benchmarks/results/actual_17_screens.json').read_text())
    assert small['status']==large['status']=='PASS'
    rows=[]
    for c in small['cases']:
        a,b=c['variants'][0],c['variants'][-1]
        rows.append({'instance':c['name'],'n':c['parameters']['n'],'k':c['parameters']['k'],
                     'target_d':c['parameters']['target_distance'],'baseline_complete_maps':a['complete_maps_tested'],
                     'filtered_complete_maps':b['complete_maps_tested'],'valid_maps':a['solutions'],
                     'baseline_median_seconds':a['median_seconds'],'filtered_median_seconds':b['median_seconds'],
                     'identical_solution_sets':c['all_variants_identical_solution_sets']})
    with (ROOT/'benchmarks/results/small_instance_comparison.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    lines=['# Measured effects of the quantum extension framework','',
           '31 October 2026. New comparisons performed on the restored project in `Desktop/quantum`. These experiments measure the existing extension framework; they do not establish a new quantum-code nonexistence result.','',
           '## Results on the actual seventeen-qubit case lists','',
           'The two q = 3 branches contain 118 and 5,995 shortened-code classes. Every physical kernel plane was regenerated from its radical. The regenerated planes exactly matched the supplied orbit partitions, and the two necessary rejection conditions were recomputed on every plane, including every member of each orbit.','',
           '| Branch | Search tasks before these tests, after symmetry | Search tasks remaining | Tasks avoided | Reduction |',
           '|---|---:|---:|---:|---:|']
    for b in large['branches']:
        label='q = 3, t = 0, hull dimension 5' if b['branch']=='q3_t0_h5' else 'q = 3, t = 2, hull dimension 3'
        lines.append(f"| {label} | {b['search_tasks_before_screens_after_symmetry']:,} | {b['search_tasks_after_screens']:,} | {b['search_tasks_avoided']:,} | {b['percent_search_tasks_avoided']:.1f}% |")
    lines += ['',
        'Here one task is a specified local type and a kernel orbit that would otherwise be passed to the extension search. There are five local types per surviving orbit in the hull-three branch and one in the hull-five branch. These task counts are not counts of complete codes or of individual search-tree nodes.','',
        'The first rejection condition uses the proved restriction that there are at most two stabilizer words of weights six or seven. A kernel plane with three such words is impossible. The second identifies support patterns already covered by the rank-two argument, so those patterns need not be searched again in the exclusive rank-one branch. This is a division of the proof cases, rather than a new independent exclusion theorem.','',
        'The reductions are 44.9% and 50.7% of the tasks that remain after symmetry identification. No runtime saving for these full branches was measured; the difficult cases need not all cost the same amount of time.','',
        '![Search tasks before and after the two necessary tests](benchmarks/results/seventeen_qubit_comparison.png)','',
        '## Complete matched small-instance experiment','',
        'Eight fixed local/shortened-code extension problems were exhaustively solved in four versions: full-map enumeration; early kernel tests; early kernel and distance tests; and those tests plus the applicable rank, hull and counting restrictions. Each version uses the same final predicate from the original mathematical model. For every instance, all four versions returned exactly the same set of accepted matrices, compared by the complete sets and their SHA-256 digests.','',
        'Every accepted map was additionally reconstructed as a physical Pauli space. The checks verify the required dimensions, commutation, exact local stabilizer and the required minimum weight outside that local stabilizer. There were 2,124 accepted-map checks across the eight instances.','',
        '| n | k | Target d | Complete maps: baseline | Complete maps: all tests | Valid maps | Baseline median (ms) | All-tests median (ms) |',
        '|---:|---:|---:|---:|---:|---:|---:|---:|']
    for row in rows:
        lines.append('| '+ ' | '.join([str(row['n']),str(row['k']),str(row['target_d']),f"{row['baseline_complete_maps']:,}",f"{row['filtered_complete_maps']:,}",str(row['valid_maps']),f"{row['baseline_median_seconds']*1000:.3f}",f"{row['filtered_median_seconds']*1000:.3f}"])+ ' |')
    lines += ['',
        'Each timing is the median of nine completed runs on the same computer. Counts are deterministic; timings are affected by system load. The shortest runs are too brief for reliable speedup ratios. The simplest feasible instance can be slower with the additional tests, because the tests themselves have a cost. Zero complete-map checks means rejection before a complete matrix was evaluated.','',
        'The k = 0 rows concern stabilizer-state extension problems. In those rows the target is the stated weight condition outside the prescribed local stabilizer; no claim about an optimal logical-code distance is made. All cases use four local sites with the balanced local stabilizer generated by XXXX and ZZZZ. The precise shortened-code generators are in the result file.','',
        'The baseline is straightforward exhaustive enumeration in Python, not the best available quantum-code algorithm. Consequently these timings do not establish superiority over another researcher\'s algorithm or predict a speedup for the large sixteen- or seventeen-qubit searches.','',
        '## Reproduce the comparisons','',
        'From the project root, after installing `requirements.txt`:','',
        '```sh',
        'python benchmarks/compare_framework.py --repeats 9',
        'python benchmarks/measure_17_screens.py',
        'python benchmarks/write_report.py',
        '```','',
        'The JSON results preserve all four variants, all repetition times, exact solution-set hashes, full instance parameters, and hashes of the input plans. The CSV is a convenient summary. The figure is a visualization of those counts.','',
        '## What remains outside these experiments','',
        'These are measurements and consistency checks of the existing framework. They do not replace the analytic proof, an external expert review, the large extension searches, or a fresh re-execution of those searches. The separate reproducibility guide records exactly which saved evidence was checked and how to run the deeper semantic replay.']
    (ROOT/'BENCHMARK_REPORT.md').write_text('\n'.join(lines)+'\n')
    print('Wrote BENCHMARK_REPORT.md and CSV summary.')
if __name__=='__main__':main()
