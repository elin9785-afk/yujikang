import random
from collections import Counter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from dna_tools import (gc, frame_stats, six_frames, metrics, shuffle,
                       merge_intervals, complement, validate_cds, spearman)
from analyze import analyze_cds, region_stats

n = 0


def check(name, actual, expected):
    global n
    if actual != expected:
        raise AssertionError(f'{name}: {actual!r} != {expected!r}')
    n += 1
    print('PASS:',name)


check('GC excludes N',gc('ATGCN'),50)
check('GC no valid base',gc('NNN'),'')
check('two stops one codon between',frame_stats('TAACCCTAG',0),dict(valid=3,stops=2,gaps=1,gap_sum=1))
check('adjacent stops gap zero',metrics(frame_stats('TAATAG',0))['mean_gap'],0)
check('no complete interval',metrics(frame_stats('ATGTAA',0))['mean_gap'],'')
check('ambiguous codon breaks gap',frame_stats('TAANNNTAG',0)['gaps'],0)
check('offset 1',frame_stats('ATAACCC',1)['stops'],1)
check('table 4 TGA not stop',frame_stats('TGATAG',0,4)['stops'],1)
check('six frames',len(six_frames('ATGAAATAA')),6)
check('reverse strand stop',six_frames('TTACAT')[3]['stops'],1)
s = 'AACGTNNGCAA'
t = shuffle(s,random.Random(42))
check('shuffle counts',Counter(t),Counter(s))
check('ambiguous positions',t[5:7],'NN')
check('shuffle deterministic',t,shuffle(s,random.Random(42)))
check('merge overlap',merge_intervals([(0,10),(5,15),(18,20)]),[[0,15],[18,20]])
check('complement',complement([[0,15],[18,20]],25),[[15,18],[20,25]])
check('normal CDS',validate_cds('ATGAAATAA'),('ok',2))
check('alternative start',validate_cds('GTGAAATAA'),('ok',2))
check('internal stop',validate_cds('ATGTAATAA')[0],'internal_stop')
check('bad length',validate_cds('ATGA')[0],'bad_length')
check('missing terminal stop',validate_cds('ATGAAA')[0],'missing_stop')
check('ambiguous CDS',validate_cds('ATGNNNTAA')[0],'empty_or_ambiguous')
check('positive rank correlation',spearman([1,2,3],[2,4,8]),1)
check('negative rank correlation',spearman([1,2,3],[8,4,2]),-1)
check('rank ties',spearman([1,1,3],[4,4,8]),1)
check('constant rank undefined',spearman([1,1],[2,3]),'')
# Two separated intervals contain just one stop each; never join them.
check('region boundaries',region_stats('TAAAAATAG',[[0,3],[6,9]],'x','CDS',11)['gaps'],0)
r = SeqRecord(Seq('TTATTTCATATGAAATAA'),id='test')
r.features = [SeqFeature(SimpleLocation(0,9,strand=-1),type='CDS',qualifiers={'locus_tag':['a']}),
              SeqFeature(SimpleLocation(9,18),type='gene',qualifiers={'locus_tag':['b'],'pseudo':['']}),
              SeqFeature(SimpleLocation(9,18),type='CDS',qualifiers={'locus_tag':['b']})]
intervals,genes,shifts = analyze_cds(r,'test',11)
check('reverse CDS',genes[0]['status'],'ok')
check('annotated pseudo excluded',intervals,[[0,9]])
check('normal reading frame internal stops',shifts[0]['stops'],0)
print(n,'checks passed')
