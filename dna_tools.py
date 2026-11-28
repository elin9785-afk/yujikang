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
                # 두 종결 코돈 사이에 있는 코돈만 센다.
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
    # A, C, G, T만 섞어서 N 같은 모호한 염기는 원래 자리에 둔다.
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
