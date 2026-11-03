"""Freeze manifest, union inputs and full-group files for currently unresolved 888 mu6 class triples."""
import json
from mixed_frames.stage4.run import read_manifest
from .common import (BINARY, BINARY_SHA, HERE, PROFILES, atomic_json, checkpoint_status, group_text,
                     legacy_rows, own_hashes, sha, union_text)


def main():
    assert sha(BINARY) == BINARY_SHA
    read_manifest()  # frozen stage-4 verification, implementation hashes, 148 inherited exclusions
    rows = legacy_rows(); (HERE / 'inputs').mkdir(exist_ok=True); classes = []
    for i in range(56):
        firsts = rows[f'888_{i}_10']['first_images']
        assert all(rows[f'888_{i}_{pf}']['first_images'] == firsts for pf in PROFILES)
        if all(checkpoint_status(f'888_{i}_{pf}', f) == 'COMPLETE' for pf in PROFILES for f in firsts):
            continue  # fully legacy now; still re-verified root-by-root by the strict aggregate
        inp = HERE / f'inputs/union_{i}.txt'; inp.write_text(union_text(i))
        text, order = group_text(i); grp = None
        if text is not None:
            grp = HERE / f'inputs/group_{i}.txt'; grp.write_text(text)
        mask = rows[f'888_{i}_10']['symmetry']
        classes.append(dict(class_index=i, first_images=firsts, input=str(inp.relative_to(HERE)), input_sha256=sha(inp),
                            group=grp and str(grp.relative_to(HERE)), group_sha256=grp and sha(grp), group_order=order,
                            mask=mask, mask_sha256=sha(HERE.parent / 'mixed_frames/stage4' / mask),
                            legacy_keys=[f'888_{i}_{pf}' for pf in PROFILES]))
    atomic_json(HERE / 'manifest.json', dict(classes=classes, hashes=own_hashes(),
        scope='888 profiles (1,0),(1,1),(2,0) superset; union budget (2,2,1); full suffix group x S8 joint-lex; '
              'Y by frozen fixed_z.solve(find_all=False), no limits.'))
    print(json.dumps(dict(classes=[c['class_index'] for c in classes], roots=sum(len(c['first_images']) for c in classes))))


if __name__ == '__main__':
    main()
