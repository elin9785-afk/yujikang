import hashlib
from pathlib import Path
from time import sleep
from urllib.parse import urlencode
from urllib.request import urlopen, Request

# Download only missing files. Never overwrite existing original files.
# EFetch may return changed annotations; reject any byte mismatch.
for line in Path("data/input_sha256.txt").read_text().splitlines():
    expected, name = line.split(maxsplit=1)
    target = Path(name)
    if target.exists():
        if hashlib.sha256(target.read_bytes()).hexdigest() != expected:
            raise SystemExit(f"Checksum mismatch: {target}")
        print("Verified:", target)
        continue
    accession = target.stem
    rettypes = ["gbwithparts", "gb"] if target.suffix == ".gb" else ["fasta"]
    for rettype in rettypes:
        query = urlencode({"db": "nuccore", "id": accession,
                           "rettype": rettype, "retmode": "text",
                           "tool": "gc_stop_tutorial"})
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + query
        with urlopen(Request(url, headers={"User-Agent": "GC-stop-tutorial/1.0"}), timeout=180) as response:
            payload = response.read()
        sleep(0.5)
        if hashlib.sha256(payload).hexdigest() == expected:
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
            print("Verified and saved:", target)
            break
    else:
        raise SystemExit(f"NCBI response differs from original: {target}. Obtain the original uploaded input.")
