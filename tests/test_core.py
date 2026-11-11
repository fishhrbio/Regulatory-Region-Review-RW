"""Mathematical reference tests independent of interval-sweep implementation."""
import random
import pytest
from chromatin_atlas.core import Interval, Atom, validate_intervals, intersect_tracks, eligible_atoms, class_overlap, summarize, read_track


def test_half_open_example():
    atoms=intersect_tracks([Interval(100,200,1)],[Interval(150,250,1)])
    assert atoms == [Atom(150,200,1,1)]
    # Full union for the hand example includes singly covered flanks.
    shared=sum(a.bp for a in atoms)
    assert shared==50
    assert shared/(100+100-shared)==pytest.approx(1/3)
    # Within COMMON territory both annotations agree completely.
    assert class_overlap(atoms,[1])["jaccard"]==1


def test_touching_and_gaps():
    assert intersect_tracks([Interval(0,10,1)],[Interval(10,20,1)])==[]
    assert intersect_tracks([Interval(0,5,1),Interval(10,15,4)],[Interval(3,12,5)])==[Atom(3,5,1,5),Atom(10,12,4,5)]


@pytest.mark.parametrize('rows',[[Interval(2,2,1)],[Interval(-1,2,1)],[Interval(0,2,16)], [Interval(0,10,1),Interval(9,12,2)], [Interval(0,10,1),Interval(0,10,1)], [Interval(10,15,1),Interval(0,5,2)]])
def test_bad_intervals_fail(rows):
    with pytest.raises(ValueError): validate_intervals(rows)


def test_mask_applies_to_both_lines_and_state13_is_valid():
    atoms=[Atom(0,10,14,1),Atom(10,20,1,15),Atom(20,25,13,4),Atom(25,28,1,1)]
    kept=eligible_atoms(atoms,[14,15])
    assert sum(a.bp for a in kept)==8
    assert class_overlap(kept,[4,5])["k_only_bp"]==5
    assert class_overlap(kept,[1])["shared_bp"]==3


def test_undefined_empty_class_and_empty_territory():
    assert class_overlap([Atom(0,10,13,13)],[1])["jaccard"] is None
    assert class_overlap([],[1])["gm_retained"] is None
    with pytest.raises(ValueError): summarize([],{'p':[1]},[])


def test_state4_and_state5_are_same_class():
    r=class_overlap([Atom(0,7,4,5),Atom(7,12,1,4),Atom(12,15,4,13)],[4,5])
    assert (r['gm_only_bp'],r['shared_bp'],r['k_only_bp'],r['union_bp'])==(3,7,5,15)
    assert r['jaccard']==pytest.approx(7/15)


def test_randomized_sweep_against_per_base_oracle():
    rng=random.Random(20261109)
    for _ in range(200):
        tracks=[]
        for line in range(2):
            x=0; rows=[]
            while x<200:
                x+=rng.randint(0,8)
                end=x+rng.randint(1,14)
                rows.append(Interval(x,end,rng.randint(1,15)))
                x=end
            tracks.append(rows)
        maps=[{p:r.state for r in rows for p in range(r.start,r.end)} for rows in tracks]
        expected={(p,maps[0][p],maps[1][p]) for p in maps[0].keys()&maps[1].keys()}
        atoms=intersect_tracks(*tracks)
        actual={(p,a.gm,a.k) for a in atoms for p in range(a.start,a.end)}
        assert actual==expected
        assert sum(a.bp for a in atoms)==len(expected)
        for excluded in [[],[14,15]]:
            kept={(p,g,k) for p,g,k in expected if g not in excluded and k not in excluded}
            if not kept: continue
            result=summarize(atoms,{'enhancer':[4,5]},excluded)
            assert result['eligible_bp']==len(kept)
            assert sum(x['gm_bp'] for x in result['coverage'])==len(kept)
            assert sum(x['k_bp'] for x in result['coverage'])==len(kept)
            r=result['overlap']['enhancer']
            shared=sum(g in [4,5] and k in [4,5] for _,g,k in kept)
            union=sum(g in [4,5] or k in [4,5] for _,g,k in kept)
            assert r['shared_bp']==shared and r['union_bp']==union
            assert r['jaccard']==(shared/union if union else None)
        reversed_atoms=intersect_tracks(tracks[1],tracks[0])
        assert class_overlap(atoms,[1])['jaccard']==class_overlap(reversed_atoms,[1])['jaccard']


def test_parser_rejects_missing_label_and_duplicates(tmp_path):
    import gzip
    path=tmp_path/'bad.bed.gz'
    for content in ['chr22\t0\t10\t\n','chr22\t0\t10\t1_Active\nchr22\t0\t10\t1_Active\n']:
        with gzip.open(path,'wt') as f: f.write(content)
        with pytest.raises(ValueError): read_track(path,'chr22')


def test_marginals_match_direct_interval_lengths():
    gm=[Interval(0,10,1),Interval(10,30,4)]
    k=[Interval(0,5,1),Interval(5,30,13)]
    r=summarize(intersect_tracks(gm,k),{'p':[1],'e':[4,5]},[])
    assert r['coverage'][0]['gm_bp']==10
    assert r['coverage'][3]['gm_bp']==20
    assert r['coverage'][12]['k_bp']==25
    assert r['overlap']['p']['jaccard']==.5


def test_pipeline_rejects_mislabeled_configuration_before_reading_data(tmp_path):
    import json
    from chromatin_atlas.pipeline import run
    config={'assembly':'hg19','chromosome':'chr22','primary_excluded_states':[14,15],
            'classes':{'active_promoter':[2],'strong_enhancer':[4,5],'broad_enhancer':[4,5,6,7]}}
    path=tmp_path/'config.json'
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError,match='Update labels'):
        run(tmp_path,path,tmp_path,False)
