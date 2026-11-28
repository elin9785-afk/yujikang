# 다음 종결 코돈은 얼마나 멀리 있을까?

**DNA의 여섯 읽기 틀에서 살펴본 GC 함량과 종결 코돈 간격**

충남대학교 생명정보융합학과 강유지

DNA 서열을 FASTA로 입력하면 양쪽 가닥을 각각 세 위치에서 읽고, 같은 읽기 틀에서 연속한 두 종결 코돈 사이에 몇 개의 코돈이 있는지 계산한다. 이 프로그램을 세균 9종의 대표 염색체에 적용해 GC 함량과 평균 간격의 관계를 조사했다.

## 내 FASTA 분석하기

<details>
<summary>처음 실행할 때 · 저장소 복제와 환경 설정</summary>

```bash
git clone https://github.com/elin9785-afk/yujikang.git
cd yujikang
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-verified.txt
```

</details>

저장소 최상위 폴더에서 `my_sequence.fasta`를 원하는 파일 경로로 바꿔 실행한다. 균주 목록에 등록할 필요가 없다.

```bash
python analyze.py --fasta-only my_sequence.fasta --out my_result
head -n 7 my_result/fasta_frames.csv
```

각 FASTA 레코드를 독립적으로 읽어 양쪽 가닥의 세 읽기 틀을 계산한다. 레코드 하나당 결과는 여섯 행이다. GC 함량, 종결 코돈 개수(`stops`), 완전한 종결 간격 수(`gaps`), 평균 간격(`mean_gap`)을 확인할 수 있다. 기본 번역표는 세균용 11번이다. CDS 안팎 비교에는 GenBank 주석이 필요하다.

손으로 검산하려면 아래 내용을 `my_sequence.fasta`로 저장해 실행한다.

```text
>demo
TAAAAGGGCTAG
```

+1 틀은 `TAA | AAG | GGC | TAG`이므로 기대값은 `stops=2`, `gaps=1`, `mean_gap=2`이다. 가운데에 `CCC`를 하나 추가한 `TAAAAGGGCCCCTAG`을 넣으면 평균 간격은 3이 된다. 계산할 수 있는 종결 간격이 없는 틀에서는 평균 칸이 비어 있다.

실제 염색체로 시연하려면 입력 경로에 저장소의 [Buchnera FASTA](data/NC_002528.1.fasta)를 지정한다.

## 코드가 계산하는 값

[analyze.py](analyze.py)가 입력을 읽고 [dna_tools.py](dna_tools.py)가 읽기 틀별 종결 개수와 간격을 센다. 전체 평균은 여섯 틀에서 얻은 모든 완전한 간격을 합쳐 계산한다. 첫 종결 앞과 마지막 종결 뒤, 모호한 코돈을 가로지르는 구간, 염색체 양 끝을 잇는 구간은 제외한다. GC 함량은 A·C·G·T의 총수를 분모로 계산한다.

| 질문 | 계산 방법 | 결과 |
|---|---|---|
| GC 함량과 간격이 함께 변하나? | 각 염색체의 여섯 읽기 틀을 조사하고 9종의 순위를 비교 | [균주별 요약](results/strain_summary.csv), [상관계수](results/correlations.csv) |
| 같은 조성에서도 서열 순서가 영향을 주나? | 각 염기의 수와 모호한 염기의 위치를 유지한 채 세 번 섞어 평균 비교 | [셔플 대조군](results/shuffle.csv) |
| CDS 안팎에서 차이가 있나? | pseudo 표시를 제외한 CDS 좌표를 합친 뒤, CDS와 나머지 구간을 각각 여섯 틀로 조사 | [영역별 빈도와 간격](results/regions.csv) |
| 주석된 단백질 길이는 어떤가? | 시작·종결과 내부 종결 유무 등을 검사한 CDS에서 아미노산 수 계산 | [CDS별 결과](results/genes.csv) |
| CDS를 한두 염기 옮겨 읽으면 어떻게 되나? | 검사를 통과한 CDS에서 말단 종결을 떼고 0, 1, 2염기 이동 후 재계산 | [읽기 틀 이동](results/shifts.csv) |

비CDS는 여기서 정한 CDS 좌표 밖의 구간이다. RNA와 조절 서열도 포함될 수 있다. 여섯 틀을 훑어 얻은 간격은 DNA 서열의 특성이며, 실제로 번역된 단백질 길이는 별도의 CDS 지표로 비교했다.

## 세균 9종에 적용한 결과

![세균 9종의 GC 함량과 평균 종결 간격](results/figures/02_gc_gap.png)

| Species | GC (%) | 평균 종결 간격 (codons) | 셔플 평균 간격 | 코딩 밀도 (%) |
|---|---:|---:|---:|---:|
| S_coelicolor | 72.12 | 84.18 | 58.70 | 87.20 |
| P_aeruginosa | 66.56 | 61.82 | 41.98 | 89.19 |
| C_crescentus | 67.21 | 53.82 | 43.53 | 89.98 |
| E_coli | 50.79 | 24.83 | 20.90 | 85.75 |
| C_difficile | 29.06 | 12.57 | 11.32 | 82.48 |
| M_tuberculosis | 65.61 | 52.28 | 39.89 | 89.47 |
| M_leprae | 57.80 | 31.70 | 27.51 | 53.17 |
| R_prowazekii | 29.00 | 11.60 | 11.31 | 73.19 |
| B_aphidicola | 26.31 | 12.22 | 10.64 | 86.38 |

GC 함량과 평균 종결 간격의 Spearman 순위상관은 0.9667이었다. 9종 모두 원래 서열의 간격이 셔플 평균보다 길었다. 평균 간격과 검증된 CDS의 평균 아미노산 길이 사이의 순위상관은 0.3667이었다. 표본은 9종이고 종 사이의 계통 관계는 보정하지 않았다.

[그림 8개](results/figures/) · [결과표 설명](results/README.md) · [입력과 실행 설정](results/run_info.json)

영역 비교 그림의 세로축은 유효 코돈 1,000개당 종결 빈도다. 평균 간격은 `regions.csv`의 `mean_gap` 열에 있다. 정상 읽기 틀의 내부 종결이 0인 것은 CDS 선별 조건에 따른 결과다.

## 전체 분석 재현하기

[data](data/)에는 8종의 FASTA·GenBank를 나누어 담은 ZIP 3개와, 시연용으로 풀어 둔 나머지 1종의 두 파일이 있다. 아래 명령은 로컬 원본을 복원하고 [18개 체크섬](data/input_sha256.txt)을 검사한 뒤 새 폴더에 결과를 쓴다.

```bash
python download_data.py
sha256sum -c data/input_sha256.txt
bash run_all.sh strains.csv data reproduced_results 3
```

[run_all.sh](run_all.sh)는 [테스트](test_all.py), [유전체 분석](analyze.py), [종 간 비교](compare.py), [그림 생성](plot.py), [출력 검사](check_output.py) 순으로 실행한다. `--fasta-only`는 서열만 읽고, CDS 분석에는 GenBank 주석을 사용한다. 다른 GenBank 입력은 `strains.csv`와 같은 열을 가진 별도 목록에 파일명과 번역표를 지정한다.

**재현성 검사** — 아래 명령은 임시 폴더에서 재실행한 결과와 저장된 표를 대조한다.

```bash
python verify_reproducibility.py
python verify_reproducibility.py --full
```

전체 검사에서 원본 18개의 체크섬, 결과 CSV 8개의 모든 행, 입력 해시와 난수 설정(seed 42, 셔플 3회)이 일치했다. Python 3.12 환경에서도 기존 Python 3.13.5 결과와 CSV 값이 같았다. [기존 실행 기록](verification/reproduce_check.txt)

<details>
<summary>초기 2종 분석에서 달라진 점</summary>

| 균주 | 초기 입력 길이 (bp) | 초기 GC / 거리 | 최종 염색체 길이 (bp) | 최종 GC / 간격 |
|---|---:|---:|---:|---:|
| *C. difficile* | 4,298,133 | 29.06% / 13.57 | 4,290,252 | 29.06% / 12.57 |
| *S. coelicolor* | 9,054,847 | 72.00% / 84.25 | 8,667,507 | 72.12% / 84.18 |

[초기 코드](initial/step1_gc_6frame.py)는 시작 위치부터 첫 종결까지도 세고 끝의 종결 코돈을 거리에 포함한다. [초기 결과표](initial/initial_summary.csv)와 최종 분석은 입력 길이도 다르다. 초기 분석에 쓰인 두 FASTA는 이 저장소에 들어 있지 않다. 최종 분석에서는 처음 계획한 RSCU 대신 셔플 대조군으로 염기 조성을 통제했다. 이 비교만으로 자연선택의 효과를 분리할 수는 없다.

</details>

## LLM 활용 기록

ChatGPT로 역상보 가닥의 읽기 틀을 구현하는 방법과 셔플에서 보존할 조건을 공부했다. 여러 균주의 입력을 처리하는 코드 작성·수정, 문서 정리, 별도 환경의 재실행 검증에도 도움을 받았다.
