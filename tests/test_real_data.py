"""Real-source conservation and sampled per-base reference calculation."""
from pathlib import Path
import hashlib
import json
import random
from chromatin_atlas.core import read_track, intersect_tracks, eligible_atoms, class_overlap

ROOT=Path(__file__).resolve().parents[1]


def test_source_hashes_and_real_sampled_base_oracle():
    manifest=json.loads((ROOT/'INPUT_MANIFEST.json').read_text())
    tracks=[]
    for f in manifest['files']:
        path=ROOT/f['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==f['sha256']
        tracks.append(read_track(path,'chr22')[0])
    atoms=intersect_tracks(*tracks)
    saved=json.loads((ROOT/'data/processed/results.json').read_text())
    assert sum(a.bp for a in atoms)==saved['common_covered_bp']
    rng=random.Random(20261109)
    # 1,000,000 sampled bp in 100 windows. Repeated windows may overlap;
    # this is an arithmetic check, not a genomic sampling inference.
    for _ in range(100):
        start=rng.randrange(16050000,51234534)
        end=start+10000
        maps=[]
        for rows in tracks:
            values={}
            for r in rows:
                if r.start>=end: break
                if r.end>start:
                    values.update({p:r.state for p in range(max(r.start,start),min(r.end,end))})
            maps.append(values)
        covered=maps[0].keys()&maps[1].keys()
        for excluded in [[],[14,15]]:
            selected={p for p in covered if maps[0][p] not in excluded and maps[1][p] not in excluded}
            # Clip the already-computed atomic intervals for the same window.
            from chromatin_atlas.core import Atom
            clipped=[Atom(max(a.start,start),min(a.end,end),a.gm,a.k)
                     for a in atoms if a.start<end and a.end>start]
            kept=eligible_atoms(clipped,excluded)
            assert sum(a.bp for a in kept)==len(selected)
            for states in [[1],[4,5],[4,5,6,7]]:
                union={p for p in selected if maps[0][p] in states or maps[1][p] in states}
                shared={p for p in selected if maps[0][p] in states and maps[1][p] in states}
                result=class_overlap(kept,states)
                assert result['union_bp']==len(union)
                assert result['shared_bp']==len(shared)
