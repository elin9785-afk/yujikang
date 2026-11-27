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
