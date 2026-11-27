"""Restore exact original inputs from ZIP archives or NCBI EFetch."""
import hashlib
from pathlib import Path
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from zipfile import ZipFile
import argparse


def digest(payload):
    return hashlib.sha256(payload).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--ncbi', action='store_true', help='Try NCBI for missing files')
    args = p.parse_args()
    data = Path('data')
    manifest = data / 'input_sha256.txt'
    missing = []
    archives = sorted(data.glob('originals_*.zip'))
    for line in manifest.read_text().splitlines():
        expected, name = line.split(maxsplit=1)
        target = Path(name)
        if target.exists():
            if digest(target.read_bytes()) != expected:
                raise SystemExit(f'Checksum mismatch (not overwritten): {target}')
            print('OK:', target)
            continue
        payload = None
        for archive in archives:
            with ZipFile(archive) as z:
                if name in z.namelist():
                    payload = z.read(name)
                    if digest(payload) != expected:
                        raise SystemExit(f'Invalid original in {archive}: {name}')
                    break
        if payload is None and args.ncbi:
            accession = target.stem
            rettypes = ['gbwithparts', 'gb'] if target.suffix == '.gb' else ['fasta']
            for rettype in rettypes:
                query = urlencode({'db': 'nuccore', 'id': accession,
                                   'rettype': rettype, 'retmode': 'text',
                                   'tool': 'gc_stop_tutorial'})
                url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?' + query
                try:
                    with urlopen(Request(url, headers={'User-Agent': 'GC-stop-tutorial/1.0'}), timeout=180) as response:
                        candidate = response.read()
                except OSError as exc:
                    print(f'NCBI request failed for {name}: {exc}')
                    continue
                finally:
                    sleep(0.5)
                if digest(candidate) == expected:
                    payload = candidate
                    break
        if payload is None:
            missing.append(name)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(payload)
        print('Restored and verified:', target)
    if missing:
        print('Missing or not byte-identical to original:')
        print('\n'.join(missing))
        raise SystemExit('Obtain the original files or originals_*.zip archives. No mismatched data saved.')
    print('PASS: all original input checksums match')


if __name__ == '__main__':
    main()
