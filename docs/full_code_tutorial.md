# 전체 코드 실습 튜토리얼

[메인 연구 페이지](../README.md)

최종 발표용 검증 코드 전체이다. 학교 서버 소스와 동일한지는 미확인이다. 파일별 코드를 복사해 같은 이름으로 저장한다.

## Step 01. 라이브러리 설정

파일: [requirements.txt](../requirements.txt)

```text
biopython>=1.81,<2
matplotlib>=3.7,<4
```

## Step 02. 9종 메타데이터

파일: [strains.csv](../strains.csv)

```csv
strain,file,accession,lifestyle,table
S_coelicolor,NC_003888.3.gb,NC_003888.3,Free_living,11
P_aeruginosa,NC_002516.2.gb,NC_002516.2,Free_living,11
C_crescentus,NC_002696.2.gb,NC_002696.2,Free_living,11
E_coli,NC_000913.3.gb,NC_000913.3,Host_associated,11
C_difficile,NC_009089.1.gb,NC_009089.1,Host_associated,11
M_tuberculosis,NC_000962.3.gb,NC_000962.3,Host_associated,11
M_leprae,NC_002677.1.gb,NC_002677.1,Obligate_intracellular,11
R_prowazekii,NC_000963.1.gb,NC_000963.1,Obligate_intracellular,11
B_aphidicola,NC_002528.1.gb,NC_002528.1,Obligate_intracellular,11
```

## Step 03. GC·6-frame·셔플·CDS 계산 함수

파일: [dna_tools.py](../dna_tools.py)

```python
import csv
import random
from math import sqrt
from Bio.Seq import Seq
from Bio.Data import CodonTable


def read_csv(path):
    with open(path, newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns):
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        w.writerows(rows)


def ratio(a, b, scale=1):
    return scale * a / b if b else ''


def gc(s):
    n = sum(s.count(b) for b in 'ACGT')
    return ratio(s.count('G') + s.count('C'), n, 100)


def frame_stats(s, offset, table=11):
    code = CodonTable.unambiguous_dna_by_id[table]
    stop_codons = set(code.stop_codons)
    valid_codons = set(code.forward_table) | stop_codons
    valid = stops = gaps = gap_sum = 0
    previous = None
    for pos in range(offset, len(s)-2, 3):
        codon = s[pos:pos+3]
        if codon not in valid_codons:
            previous = None
            continue
        valid += 1
        if codon in stop_codons:
            stops += 1
            if previous is not None:
                gap_sum += (pos-previous)//3 - 1
                gaps += 1
            previous = pos
    return dict(valid=valid, stops=stops, gaps=gaps, gap_sum=gap_sum)


def six_frames(s, table=11):
    rows = []
    for sign, seq in [('+', s), ('-', str(Seq(s).reverse_complement()))]:
        for offset in range(3):
            row = frame_stats(seq, offset, table)
            row['frame'] = sign + str(offset+1)
            rows.append(row)
    return rows


def total_stats(rows):
    total = dict(valid=0, stops=0, gaps=0, gap_sum=0)
    for row in rows:
        for k in total:
            total[k] += row[k]
    return total


def metrics(row):
    result = dict(row)
    result['stop_rate'] = ratio(row['stops'], row['valid'], 1000)
    result['mean_gap'] = ratio(row['gap_sum'], row['gaps'])
    return result


def shuffle(s, rng):
    # Preserve A/C/G/T counts and positions of ambiguous bases.
    bases = [b for b in s if b in 'ACGT']
    rng.shuffle(bases)
    it = iter(bases)
    return ''.join(next(it) if b in 'ACGT' else b for b in s)


def merge_intervals(intervals):
    merged = []
    for a, b in sorted(intervals):
        if a >= b:
            continue
        if merged and a <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], b)
        else:
            merged.append([a, b])
    return merged


def complement(intervals, length):
    result = []
    start = 0
    for a, b in intervals:
        if a > start:
            result.append([start, a])
        start = b
    if start < length:
        result.append([start, length])
    return result


def validate_cds(s, table=11):
    if not s or any(b not in 'ACGT' for b in s):
        return 'empty_or_ambiguous', ''
    if len(s) % 3:
        return 'bad_length', ''
    protein = str(Seq(s).translate(table=table, to_stop=False))
    if '*' in protein[:-1]:
        return 'internal_stop', ''
    if not protein.endswith('*'):
        return 'missing_stop', ''
    code = CodonTable.unambiguous_dna_by_id[table]
    if len(s) < 6 or s[:3] not in code.start_codons:
        return 'invalid_start', ''
    return 'ok', len(Seq(s).translate(table=table, to_stop=True))


def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        for k in range(i, j):
            result[order[k]] = (i+1+j)/2
        i = j
    return result


def spearman(x, y):
    if len(x) != len(y) or len(x) < 2:
        return ''
    a, b = ranks(x), ranks(y)
    center = (len(a)+1)/2
    top = sum((u-center)*(v-center) for u,v in zip(a,b))
    bottom = sqrt(sum((u-center)**2 for u in a)*sum((v-center)**2 for v in b))
    return top/bottom if bottom else ''
```

## Step 04. GenBank 및 FASTA 입력과 분석

파일: [analyze.py](../analyze.py)

```python
import argparse
import hashlib
import json
import platform
import random
from pathlib import Path
from statistics import mean
from Bio import SeqIO, __version__ as bio_version
from Bio.SeqFeature import ExactPosition
from dna_tools import (read_csv, write_csv, gc, ratio, frame_stats, six_frames,
                       total_stats, metrics, shuffle, merge_intervals,
                       complement, validate_cds)

STAT = ['valid','stops','gaps','gap_sum','stop_rate','mean_gap']
SUMMARY = ['strain','accession','lifestyle','genome_bp','gc_percent','mean_gap',
           'stop_rate','shuffle_mean_gap','shuffle_mean_rate','gap_difference',
           'coding_bp','coding_density','cds_count','valid_cds','mean_aa_length']


def read_genome(path, accession):
    records = []
    for r in SeqIO.parse(path, 'genbank'):
        if not accession or r.id == accession:
            records.append(r)
    if len(records) != 1:
        raise ValueError('Expected one matching chromosome: '+str(path))
    r = records[0]
    if not len(r.seq) or not any(f.type == 'CDS' for f in r.features):
        raise ValueError('Missing sequence or CDS annotation: '+r.id)
    return r


def locus(f):
    return f.qualifiers.get('locus_tag', [str(f.location)])[0]


def pseudo(f):
    return f.type == 'pseudogene' or 'pseudo' in f.qualifiers or 'pseudogene' in f.qualifiers


def analyze_cds(r, name, table):
    pseudo_ids = set()
    for f in r.features:
        if pseudo(f):
            pseudo_ids.add(locus(f))
    intervals, genes = [], []
    shift_rows = [[], [], []]
    for f in r.features:
        if f.type != 'CDS':
            continue
        q = f.qualifiers
        aa = ''
        if f.location is None:
            raise ValueError('CDS lacks location: '+r.id)
        if any(p.ref or p.ref_db for p in f.location.parts):
            raise ValueError('Remote CDS unsupported: '+r.id)
        if locus(f) in pseudo_ids:
            status = 'annotated_pseudo'
        else:
            partial = False
            for p in f.location.parts:
                a, b = int(p.start), int(p.end)
                if not 0 <= a <= b <= len(r.seq):
                    raise ValueError('Invalid CDS coordinate: '+r.id)
                intervals.append((a,b))
                if not isinstance(p.start, ExactPosition) or not isinstance(p.end, ExactPosition):
                    partial = True
            t = int(q.get('transl_table', [table])[0])
            if (partial or q.get('codon_start',['1'])[0] != '1'
                    or any(k in q for k in ['exception','transl_except','ribosomal_slippage'])):
                status = 'partial_or_exception'
            else:
                s = str(f.extract(r.seq)).upper()
                status, aa = validate_cds(s,t)
                if status == 'ok':
                    for offset in range(3):
                        shift_rows[offset].append(frame_stats(s[:-3],offset,t))
        genes.append(dict(strain=name,locus=locus(f),length_bp=len(f),status=status,aa_length=aa))
    shifts = []
    for i in range(3):
        row = metrics(total_stats(shift_rows[i]))
        row.update(strain=name,offset=i)
        shifts.append(row)
    return merge_intervals(intervals), genes, shifts


def region_stats(seq, intervals, name, region, table):
    counts = []
    bases = dict(A=0,C=0,G=0,T=0)
    bp = 0
    for a,b in intervals:
        part = seq[a:b]
        bp += len(part)
        for base in bases:
            bases[base] += part.count(base)
        # No codon or gap crosses an interval boundary.
        counts.extend(six_frames(part,table))
    result = metrics(total_stats(counts))
    result.update(strain=name,region=region,bp=bp,segments=len(intervals),
                  gc_percent=ratio(bases['G']+bases['C'],sum(bases.values()),100))
    return result


def sequence_analysis(seq, name, table, reps, seed):
    observed = None
    frames, controls = [], []
    rng = random.Random(seed)
    for rep in range(reps+1):
        s = seq if rep == 0 else shuffle(seq,rng)
        raw = six_frames(s,table)
        for f in raw:
            row = metrics(f)
            row.update(strain=name,rep=rep)
            frames.append(row)
        if rep == 0:
            observed = metrics(total_stats(raw))
        else:
            row = metrics(total_stats(raw))
            row.update(strain=name,rep=rep)
            controls.append(row)
    gaps = [r['mean_gap'] for r in controls if r['mean_gap'] != '']
    rates = [r['stop_rate'] for r in controls if r['stop_rate'] != '']
    null_gap = mean(gaps) if gaps else ''
    null_rate = mean(rates) if rates else ''
    diff = observed['mean_gap']-null_gap if observed['mean_gap'] != '' and null_gap != '' else ''
    summary = dict(strain=name,genome_bp=len(seq),gc_percent=gc(seq),
                   mean_gap=observed['mean_gap'],stop_rate=observed['stop_rate'],
                   shuffle_mean_gap=null_gap,shuffle_mean_rate=null_rate,gap_difference=diff)
    return summary, frames, controls


def main(meta,data,out,reps,seed):
    out.mkdir(parents=True,exist_ok=True)
    summaries, genes, shifts, regions, frames, controls, inputs = [],[],[],[],[],[],[]
    seen = set()
    for row in read_csv(meta):
        name = row['strain']
        if name in seen:
            raise ValueError('Duplicate strain: '+name)
        seen.add(name)
        path = data/row['file']
        r = read_genome(path,row['accession'])
        seq = str(r.seq).upper()
        table = int(row['table'])
        print('Analyzing:',name,flush=True)
        s,fr,ct = sequence_analysis(seq,name,table,reps,seed)
        intervals,g,sh = analyze_cds(r,name,table)
        reg = [region_stats(seq,intervals,name,'CDS',table),
               region_stats(seq,complement(intervals,len(seq)),name,'non_CDS',table)]
        coding_bp = sum(b-a for a,b in intervals)
        lengths = [x['aa_length'] for x in g if x['status']=='ok']
        s.update(accession=r.id,lifestyle=row['lifestyle'],coding_bp=coding_bp,
                 coding_density=ratio(coding_bp,len(seq),100),cds_count=len(g),
                 valid_cds=len(lengths),mean_aa_length=mean(lengths) if lengths else '')
        summaries.append(s)
        genes.extend(g); shifts.extend(sh); regions.extend(reg)
        frames.extend(fr); controls.extend(ct)
        inputs.append(dict(strain=name,accession=r.id,file=row['file'],table=table,
                           sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
    if not summaries:
        raise ValueError('Empty metadata')
    write_csv(out/'strain_summary.csv',summaries,SUMMARY)
    write_csv(out/'genes.csv',genes,['strain','locus','length_bp','status','aa_length'])
    write_csv(out/'frames.csv',frames,['strain','rep','frame']+STAT)
    write_csv(out/'shuffle.csv',controls,['strain','rep']+STAT)
    write_csv(out/'regions.csv',regions,['strain','region','bp','segments','gc_percent']+STAT)
    write_csv(out/'shifts.csv',shifts,['strain','offset']+STAT)
    info = dict(seed=seed,reps=reps,python=platform.python_version(),biopython=bio_version,inputs=inputs)
    (out/'run_info.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf-8')


def fasta_only(path,out,table):
    out.mkdir(parents=True,exist_ok=True)
    rows = []
    for r in SeqIO.parse(path,'fasta'):
        s = str(r.seq).upper()
        if not s or any(b not in 'ACGTRYSWKMBDHVN' for b in s):
            raise ValueError('Invalid DNA FASTA: '+r.id)
        for f in six_frames(s,table):
            row = metrics(f)
            row.update(record=r.id,genome_bp=len(s),gc_percent=gc(s))
            rows.append(row)
    if not rows:
        raise ValueError('Empty FASTA')
    write_csv(out/'fasta_frames.csv',rows,['record','genome_bp','gc_percent','frame']+STAT)
    print('Saved',out/'fasta_frames.csv')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--metadata',default='strains.csv')
    p.add_argument('--data',default='data')
    p.add_argument('--out',default='results')
    p.add_argument('--reps',type=int,default=3)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--fasta-only')
    p.add_argument('--table',type=int,default=11)
    a = p.parse_args()
    if a.fasta_only:
        fasta_only(Path(a.fasta_only),Path(a.out),a.table)
    else:
        if a.reps < 1:
            p.error('--reps must be >= 1')
        main(Path(a.metadata),Path(a.data),Path(a.out),a.reps,a.seed)
```

## Step 05. 균주 비교 및 순위상관

파일: [compare.py](../compare.py)

```python
import sys
from pathlib import Path
from statistics import mean
from dna_tools import read_csv, write_csv, spearman


def main(out):
    rows = read_csv(out/'strain_summary.csv')
    pairs = [('gc_percent','mean_gap'),('gc_percent','stop_rate'),
             ('mean_gap','mean_aa_length'),('genome_bp','coding_density')]
    correlations = []
    for a,b in pairs:
        x,y = [],[]
        for r in rows:
            if r[a] != '' and r[b] != '':
                x.append(float(r[a])); y.append(float(r[b]))
        correlations.append(dict(x=a,y=b,n=len(x),rho=spearman(x,y)))
    write_csv(out/'correlations.csv',correlations,['x','y','n','rho'])
    groups = {}
    for r in rows:
        groups.setdefault(r['lifestyle'],[]).append(r)
    result = []
    for group in groups:
        for metric in ['genome_bp','gc_percent','mean_gap','coding_density','mean_aa_length']:
            values = [float(r[metric]) for r in groups[group] if r[metric] != '']
            result.append(dict(lifestyle=group,metric=metric,n=len(values),mean=mean(values) if values else ''))
    write_csv(out/'lifestyle.csv',result,['lifestyle','metric','n','mean'])


if __name__ == '__main__':
    main(Path(sys.argv[1]))
```

## Step 06. 그래프 생성

파일: [plot.py](../plot.py)

```python
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from dna_tools import read_csv

out = Path(sys.argv[1])
fig = out/'figures'
fig.mkdir(exist_ok=True)
rows = sorted(read_csv(out/'strain_summary.csv'),key=lambda r: float(r['gc_percent']) if r['gc_percent'] else -1)


def save(name):
    plt.tight_layout()
    plt.savefig(fig/(name+'.png'),dpi=160)
    plt.close()


def bars(names,series,ylabel,name):
    plt.figure(figsize=(10,5))
    width = .8/len(series)
    for j,(label,values) in enumerate(series):
        x = []; y = []
        for i,v in enumerate(values):
            if v != '':
                x.append(i+(j-(len(series)-1)/2)*width)
                y.append(float(v))
        plt.bar(x,y,width=width,label=label)
    plt.xticks(range(len(names)),names,rotation=45,ha='right')
    plt.ylabel(ylabel)
    if len(series)>1:
        plt.legend()
    save(name)


def scatter(xkey,ykey,xlabel,ylabel,name,xscale=1):
    plt.figure(figsize=(9,6))
    groups = {}
    for r in rows:
        if r[xkey] != '' and r[ykey] != '':
            groups.setdefault(r['lifestyle'],[]).append(r)
    for group in groups:
        x = []; y = []
        for r in groups[group]:
            x.append(float(r[xkey])/xscale); y.append(float(r[ykey]))
            plt.annotate(r['strain'],(x[-1],y[-1]),fontsize=7,xytext=(3,3),textcoords='offset points')
        plt.scatter(x,y,label=group)
    plt.xlabel(xlabel); plt.ylabel(ylabel)
    if groups:
        plt.legend(fontsize=8)
    save(name)


def series_from_file(filename,key,labels,metric):
    data = read_csv(out/filename)
    series = []
    for label in labels:
        values = []
        for r in rows:
            value = ''
            for d in data:
                if d['strain']==r['strain'] and d[key]==label and d.get('rep','0')=='0':
                    value = d[metric]
                    break
            values.append(value)
        series.append((label,values))
    return series


names = [r['strain'] for r in rows]
bars(names,[('GC',[r['gc_percent'] for r in rows])],'GC (%)','01_gc')
scatter('gc_percent','mean_gap','GC (%)','Mean stop interval (codons)','02_gc_gap')
bars(names,[('Observed',[r['mean_gap'] for r in rows]),('Shuffled mean',[r['shuffle_mean_gap'] for r in rows])],
     'Mean stop interval (codons)','03_shuffle')
bars(names,series_from_file('frames.csv','frame',['+1','+2','+3','-1','-2','-3'],'mean_gap'),
     'Mean stop interval (codons)','04_frames')
bars(names,series_from_file('regions.csv','region',['CDS','non_CDS'],'stop_rate'),
     'Stops per 1000 valid codons (six frames)','05_regions')
bars(names,series_from_file('shifts.csv','offset',['0','1','2'],'stop_rate'),
     'Stops per 1000 valid codons','06_shifts')
plt.figure(figsize=(10,5))
genes = read_csv(out/'genes.csv')
values, labels = [],[]
for r in rows:
    lengths = [int(g['aa_length']) for g in genes if g['strain']==r['strain'] and g['status']=='ok']
    if lengths:
        values.append(lengths); labels.append(r['strain'])
if values:
    plt.boxplot(values,showfliers=False)
    plt.xticks(range(1,len(labels)+1),labels,rotation=45,ha='right')
else:
    plt.text(.5,.5,'No validated CDS',ha='center',transform=plt.gca().transAxes)
plt.ylabel('Validated CDS length (amino acids; outliers hidden)')
save('07_cds_lengths')
scatter('genome_bp','coding_density','Chromosome size (Mb)','Non-pseudo CDS coverage (%)','08_size_density',1e6)
```

## Step 07. 29개 검사

파일: [test_all.py](../test_all.py)

```python
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
```

## Step 08. 합성 유전체 생성

파일: [demo.py](../demo.py)

```python
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation
from dna_tools import write_csv

p = Path('demo_data')
p.mkdir(exist_ok=True)
rows = []
for i, middle in enumerate(['ATGAAATAA','ATGGCCTAA','ATGGGCGCCTAA']):
    name = 'demo'+str(i)
    # Two genuine CDS, an internal-stop CDS, and non-CDS sequence.
    first = 'ATGAAATAA'
    reverse = str(Seq(middle).reverse_complement())
    broken = 'ATGTAATAA'
    tail = 'TAACCCTAG' + 'ACGT'*(i+2)
    s = first+reverse+broken+tail
    a = len(first)
    b = a+len(reverse)
    c = b+len(broken)
    r = SeqRecord(Seq(s),id=name,description='Synthetic test only')
    r.annotations['molecule_type'] = 'DNA'
    r.annotations['date'] = '01-JAN-2000'
    r.features = [
        SeqFeature(SimpleLocation(0,a,strand=1),type='CDS',qualifiers={'locus_tag':['a']}),
        SeqFeature(SimpleLocation(a,b,strand=-1),type='CDS',qualifiers={'locus_tag':['b']}),
        SeqFeature(SimpleLocation(b,c,strand=1),type='CDS',qualifiers={'locus_tag':['c']})]
    SeqIO.write(r,p/(name+'.gb'),'genbank')
    SeqIO.write(r,p/(name+'.fasta'),'fasta')
    rows.append(dict(strain=name,file=name+'.gb',accession=name,lifestyle='Synthetic',table=11))
write_csv('demo_strains.csv',rows,['strain','file','accession','lifestyle','table'])
```

## Step 09. 출력 검증

파일: [check_output.py](../check_output.py)

```python
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
```

## Step 10. 전체 분석 실행

파일: [run_all.sh](../run_all.sh)

```bash
#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
META="${1:-strains.csv}"
DATA="${2:-data}"
OUT="${3:-results}"
REPS="${4:-3}"
python test_all.py
python analyze.py --metadata "$META" --data "$DATA" --out "$OUT" --reps "$REPS" --seed 42
python compare.py "$OUT"
python plot.py "$OUT"
python check_output.py "$META" "$OUT" "$REPS"
echo "Pipeline completed: $OUT"
```
