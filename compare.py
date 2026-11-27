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
