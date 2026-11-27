"""Rerun the published pipeline and check its newly generated results."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
SAMPLE = 'B_aphidicola'
CSV_NAMES = ('strain_summary.csv', 'frames.csv', 'genes.csv', 'regions.csv',
             'shifts.csv', 'shuffle.csv', 'correlations.csv', 'lifestyle.csv')


def require(condition, message):
    if not condition:
        raise SystemExit('FAIL: ' + message)


def check_hashes(manifest, selected=None):
    count = 0
    for line in (ROOT / manifest).read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        if selected is not None and name not in selected:
            continue
        path = ROOT / name
        require(path.is_file(), 'Missing file: ' + name)
        require(hashlib.sha256(path.read_bytes()).hexdigest() == expected,
                'Checksum mismatch: ' + name)
        count += 1
    if selected is not None:
        require(count == len(selected), 'Incomplete checksum manifest')
    print(f'PASS: {count} file checksums ({manifest})', flush=True)


def run(work, *args):
    print('RUN:', ' '.join(str(x) for x in args), flush=True)
    subprocess.run([sys.executable, *map(str, args)], cwd=work, check=True)


def rows(path):
    with path.open(newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


def pipeline(work, meta, data, out):
    run(work, ROOT / 'analyze.py', '--metadata', meta, '--data', data,
        '--out', out, '--reps', '3', '--seed', '42')
    run(work, ROOT / 'compare.py', out)
    run(work, ROOT / 'plot.py', out)
    run(work, ROOT / 'check_output.py', meta, out, '3')


def compare_csv(actual, reference, sample=None):
    for name in CSV_NAMES:
        if sample is not None and name in ('correlations.csv', 'lifestyle.csv'):
            continue
        expected = rows(reference / name)
        if sample is not None:
            expected = [row for row in expected if row['strain'] == sample]
        require(rows(actual / name) == expected, 'Different CSV: ' + name)
        print('MATCH:', name, flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--full', action='store_true', help='Restore and rerun all nine original genomes')
    args = p.parse_args()
    check_hashes('verification/source_sha256.txt')
    if args.full:
        run(ROOT, ROOT / 'download_data.py')
        check_hashes('data/input_sha256.txt')
    else:
        check_hashes('data/input_sha256.txt',
                     {'data/NC_002528.1.gb', 'data/NC_002528.1.fasta'})
    with tempfile.TemporaryDirectory(prefix='gc_reproduce_') as tmp:
        work = Path(tmp)
        run(work, ROOT / 'test_all.py')
        if not args.full:
            run(work, ROOT / 'demo.py')
            pipeline(work, 'demo_strains.csv', 'demo_data', 'demo_a')
            pipeline(work, 'demo_strains.csv', 'demo_data', 'demo_b')
            a, b = work / 'demo_a', work / 'demo_b'
            files = sorted(p.relative_to(a) for p in a.rglob('*') if p.is_file())
            require(files == sorted(p.relative_to(b) for p in b.rglob('*') if p.is_file()),
                    'Synthetic output file list differs')
            for name in files:
                require((a/name).read_bytes() == (b/name).read_bytes(),
                        'Synthetic repeat differs: ' + str(name))
            print(f'PASS: synthetic repeat ({len(files)} files)', flush=True)
            meta = work / 'single.csv'
            selected = [row for row in rows(ROOT / 'strains.csv') if row['strain'] == SAMPLE]
            require(len(selected) == 1, 'Sample metadata missing or duplicated')
            with meta.open('w', newline='', encoding='utf-8') as f:
                writer = csv.DictWriter(f, fieldnames=list(selected[0]))
                writer.writeheader()
                writer.writerows(selected)
            run(work, ROOT / 'analyze.py', '--fasta-only', ROOT / 'data/NC_002528.1.fasta',
                '--out', 'live')
            frames = rows(work / 'live/fasta_frames.csv')
            expected = [r for r in rows(ROOT / 'results/frames.csv')
                        if r['strain'] == SAMPLE and r['rep'] == '0']
            require(len(frames) == len(expected) == 6, 'Expected six frames')
            for actual, ref in zip(frames, expected):
                require(all(actual[k] == ref[k] for k in ref if k not in ('strain', 'rep')),
                        'FASTA frame differs: ' + ref['frame'])
            print('PASS: real FASTA, all six frames match', flush=True)
        else:
            meta = ROOT / 'strains.csv'
        pipeline(work, meta, ROOT / 'data', 'actual')
        compare_csv(work / 'actual', ROOT / 'results', None if args.full else SAMPLE)
        if args.full:
            actual = json.loads((work / 'actual/run_info.json').read_text())
            expected = json.loads((ROOT / 'results/run_info.json').read_text())
            for key in ('seed', 'reps', 'inputs'):
                require(actual[key] == expected[key], 'Different run setting: ' + key)
            print('PASS: original input hashes, seed and repeat count match', flush=True)
            if actual == expected:
                print('MATCH: run_info.json (including versions)', flush=True)
            else:
                print('NOTE: environment version differs; CSV values matched exactly', flush=True)
        print('PASS: FULL nine-genome reproduction' if args.full else
              'PASS: QUICK public-data reproduction (one genome + synthetic repeats)', flush=True)


if __name__ == '__main__':
    main()
