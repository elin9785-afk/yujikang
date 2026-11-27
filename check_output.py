import sys
import math
from pathlib import Path
from dna_tools import read_csv

meta, out, reps = Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3])
rows = read_csv(out/'strain_summary.csv')
assert [r['strain'] for r in rows] == [r['strain'] for r in read_csv(meta)]
assert rows
frames = read_csv(out/'frames.csv')
regions = read_csv(out/'regions.csv')
shifts = read_csv(out/'shifts.csv')
genes = read_csv(out/'genes.csv')
controls = read_csv(out/'shuffle.csv')
for r in rows:
    name = r['strain']
    f = [x for x in frames if x['strain']==name]
    assert len(f)==6*(reps+1)
    assert len({(x['rep'],x['frame']) for x in f})==len(f)
    assert len([x for x in controls if x['strain']==name])==reps
    reg = [x for x in regions if x['strain']==name]
    assert len(reg)==2
    assert sum(int(x['bp']) for x in reg)==int(r['genome_bp'])
    assert int(reg[0]['bp'])==int(r['coding_bp'])
    assert math.isclose(float(r['coding_density']),100*int(r['coding_bp'])/int(r['genome_bp']))
    assert 0 <= float(r['coding_density']) <= 100
    g = [x for x in genes if x['strain']==name]
    assert len(g)==int(r['cds_count'])
    assert len([x for x in g if x['status']=='ok'])==int(r['valid_cds'])
    assert len([x for x in shifts if x['strain']==name])==3
    obs = [x for x in f if x['rep']=='0']
    gaps = sum(int(x['gaps']) for x in obs)
    if gaps:
        assert math.isclose(float(r['mean_gap']),sum(int(x['gap_sum']) for x in obs)/gaps)
    else:
        assert r['mean_gap']==''
for filename in ['frames.csv','regions.csv','shifts.csv','shuffle.csv']:
    for r in read_csv(out/filename):
        assert 0<=int(r['gaps'])<=int(r['stops'])<=int(r['valid'])
        assert int(r['gap_sum'])>=0
assert len(read_csv(out/'correlations.csv'))==4
figs = ['01_gc','02_gc_gap','03_shuffle','04_frames','05_regions','06_shifts','07_cds_lengths','08_size_density']
for name in figs:
    assert (out/'figures'/(name+'.png')).stat().st_size>0
print('PASS: metadata order, CSV relationships, 8 figures')
