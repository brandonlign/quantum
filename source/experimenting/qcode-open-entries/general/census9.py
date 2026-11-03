import json, sys
src = open('/Users/Brandon/Desktop/research/experimenting/qcode-open-entries/release16/census/mass_census.py').read().split('if __name__')[0]
exec(src)
k = int(sys.argv[1])
found = census(9, k)
json.dump([g for g, a in found], open(f'cells9_{k}.json', 'w'))
print('saved', len(found))
