# 실제 9종 입력 데이터

학생이 첨부한 GenBank 9개 및 FASTA 9개의 accession과 염기서열 일치 여부를 확인하였다. [메타데이터](../strains.csv)와 [원본 파일 SHA-256](input_sha256.txt)을 제공한다. 원본 18개는 아직 웹 업로드가 필요하며, 업로드 시 이 폴더에 배치한다. GenBank 9개는 전체 분석에 필요하며 FASTA는 별도 시연용이다.

| accession | GenBank | FASTA |
|---|---|---|
| [NC_003888.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_003888.3) | NC_003888.3.gb | NC_003888.3.fasta |
| [NC_002516.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002516.2) | NC_002516.2.gb | NC_002516.2.fasta |
| [NC_002696.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_002696.2) | NC_002696.2.gb | NC_002696.2.fasta |
| [NC_000913.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3) | NC_000913.3.gb | NC_000913.3.fasta |
| [NC_009089.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_009089.1) | NC_009089.1.gb | NC_009089.1.fasta |
| [NC_000962.3](https://www.ncbi.nlm.nih.gov/nuccore/NC_000962.3) | NC_000962.3.gb | NC_000962.3.fasta |
| [NC_002677.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002677.1) | NC_002677.1.gb | NC_002677.1.fasta |
| [NC_000963.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_000963.1) | NC_000963.1.gb | NC_000963.1.fasta |
| [NC_002528.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_002528.1) | NC_002528.1.gb | NC_002528.1.fasta |

원본 확보 후 저장소 최상위에서 다음 명령어로 검증한다.

```bash
sha256sum -c data/input_sha256.txt
```

NCBI에서 다시 내려받는 자료는 주석 변경으로 기존 파일과 체크섬이 달라질 수 있다.
