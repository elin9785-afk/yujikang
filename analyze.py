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
