# GC content and stop codon distribution in bacterial genomes

세균의 GC 함량과 코딩 영역 특성에 따른 종결 코돈 분포 비교 분석

## Research questions
1. Can base composition alone explain stop-to-stop intervals? Compare nine chromosomes with nucleotide-preserving shuffled controls.
2. How do CDS coverage and the remaining regions differ in stop frequencies? Compare six-frame region scans and separately simulate offsets in validated CDS.
3. How do chromosome size, annotated CDS coverage and validated CDS length vary across strains?

## Inputs
`strains.csv` lists nine selected RefSeq chromosome accession versions. Download the corresponding GenBank files from https://www.ncbi.nlm.nih.gov/nuccore/ and save them in `data/` using the exact listed filenames. GenBank includes the sequence and annotation, so no additional FASTA is required. Each row analyzes one chromosome; plasmids and other records are not included. Lifestyle groups are operational sample-selection categories, not causal predictors. The same workflow supports any single matching annotated bacterial chromosome with appropriate genetic code and standard local CDS coordinates. Empty accession is allowed only when the file contains exactly one record.

## Setup (Bash / Python 3.9+)
```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m pip freeze > requirements-lock.txt
python --version > python-version.txt
```
On the original machine, an existing compatible virtual environment may be reused. For reproducibility use the same Python and package versions. Requirements ranges allow installation with compatible packages; the lock file records the actual environment.

## Run
```bash
python demo.py
bash run_all.sh demo_strains.csv demo_data demo_results 3
bash run_all.sh demo_strains.csv demo_data demo_repeat 3
diff -r demo_results demo_repeat && echo 'Reproducibility PASS'
bash run_all.sh strains.csv data pilot_results 1
bash run_all.sh strains.csv data results 3
# Optional: more shuffled controls (descriptive comparison)
bash run_all.sh strains.csv data results_20 20
```
Full chromosomes may require substantial CPU time for repeated pure Python scans. Pilot uses one shuffle, default final analysis three, 20 optional. All output paths should be distinct or empty for a clean run.

### FASTA-only demonstration
```bash
python analyze.py --fasta-only demo_data/demo0.fasta --out live
head -n 3 live/fasta_frames.csv
```
A new DNA FASTA can replace the demo path. This mode supports multiple DNA records, analyzes each separately, and reports six frames only (no CDS annotation, coding density or shuffled controls). `--table` selects the genetic code (default 11).

## Code files
- `dna_tools.py`: base composition, stop scan, shuffle, intervals, CDS validation, Spearman.
- `analyze.py`: input, six frames, control replicates, GenBank features, region and shifted-CDS analysis, CSV and input hash manifest.
- `compare.py`: exploratory Spearman rank correlation and lifestyle means.
- `plot.py`: eight figures.
- `test_all.py`: 29 hand-checkable function and feature tests.
- `demo.py`: generates three synthetic genomes (not biological data).
- `check_output.py`: checks CSV relationships, metadata order and eight nonempty figures.
- `run_all.sh`: connects tests, analysis, comparison, visualization and output checks.

## Definitions and interpretation
- GC (%) = (G+C)/(A+C+G+T)*100. Non-ACGT bases are excluded from GC and invalid triplets from stop-rate denominators.
- Each of six linear reading frames uses the specified chromosome genetic code (11 for the nine targets). A stop interval is the number of intervening codons between two successive stops in the same frame, excluding both stops. Adjacent stops give zero. Incomplete intervals at ends or across an ambiguous codon are excluded. Circular wrap-around is not analyzed. The whole-chromosome mean pools complete intervals across the six frames; it is not the mean distance from every possible start to the next stop.
- Stop rate = stops per 1000 valid triplets. Six-frame measurements are descriptive scans, not gene predictions. All six frames overlap and are not independent samples.
- Shuffle preserves the counts of A/C/G/T and ambiguous positions, but not dinucleotide frequencies, codon usage, or longer sequence features. The seed is reset per chromosome; replicate means are averaged, excluding undefined values. This controls mononucleotide composition, not GC alone in isolation or all possible confounders. Differences are descriptive, not evidence of selection by themselves.
- Coding density = union length of annotated non-pseudo CDS coordinates / selected chromosome length. Includes partial and otherwise invalid non-pseudo CDS coverage as annotated. The `non_CDS` class is the complement and may contain RNAs, pseudogenes and regulatory regions; it is not synonymous with intergenic or nonfunctional DNA. Pseudo annotations are excluded from coding coverage, not analyzed as a separate research metric. Locus tags are used to match pseudo gene/CDS annotations (location fallback when unavailable). Annotation differences can affect results.
- Region scans pool all six local frames of each merged continuous interval without joining disjoint regions. Codons and intervals cannot cross region boundaries. These region rates do not represent translated reading frames. Region GC is reported alongside rates; different GC and segment lengths can affect comparisons.
- True CDS validation uses feature extraction with strand and joined locations, feature translation table if provided, and excludes pseudo, partial and exceptional CDS from validated lengths. Standard bacterial alternative start codons are accepted. Protein length excludes the terminal stop. Offset analysis uses only validated CDS and removes the terminal stop before scanning offsets 0,1,2. Offset 0 has zero internal stops by validation design; this is not measured frameshift frequency. CDS length is not the same as six-frame stop interval.
- Spearman correlations are exploratory, with average tied ranks; blanks mean undefined, not zero. n is the number of usable chromosomes. No significance tests or phylogenetic adjustment are performed. The nine targets are not a random sample of all bacteria.
- Unsupported cases (multiple records without accession, remote CDS, missing annotations) fail explicitly. This is not a universal genome annotation tool.

## Outputs
`strain_summary.csv`, `genes.csv`, `frames.csv`, `shuffle.csv`, `regions.csv`, `shifts.csv`, `correlations.csv`, `lifestyle.csv`, `run_info.json` and `figures/`.
Eight figures: 01 GC, 02 GC vs interval, 03 observed/shuffle, 04 six frames, 05 CDS vs non-CDS stop rate, 06 offset stop rate, 07 validated CDS amino-acid length distribution, 08 chromosome size vs coding density. Blank numeric cells denote undefined metrics. Boxplot hides outlier points for legibility, not from the underlying CSV.

## Verification
Synthetic inputs and FASTA-only demonstration were checked in the development environment. Real NCBI chromosome analysis must be run on the user's files. `run_info.json` records the SHA-256 input fingerprints, genetic codes, seed, replicate count, Python and Biopython versions. `diff -r` compares outputs from repeated runs within the same environment. PNG bytes need not match across different software environments.

Development test environment: Python 3.13.5, Biopython 1.86, Matplotlib 3.10.8.
