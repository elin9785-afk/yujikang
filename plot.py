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
