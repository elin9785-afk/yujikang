# 작은 서열로 실습하기

이 폴더의 두 FASTA는 **새로 만든 검증용 합성 서열**이며 실제 세균 유전체나 학생의 기존 실험 자료가 아니다. 초기 코드가 파일명과 종 이름을 고정하므로 호환을 위해 같은 이름을 사용한다.

저장소 최상위에서 Biopython 설치 후 다음 명령어를 실행한다.

```bash
cd examples
python ../step1_gc_6frame.py
cat initial_summary.csv
```

결과는 examples/initial_summary.csv에 생성되므로 최상위의 기존 결과는 보존된다.

| 입력 | 길이(bp) | GC(%) | 평균 종결 거리(코돈) |
|---|---:|---:|---:|
| 합성 저GC 서열 | 15 | 20.0 | 1.67 |
| 합성 고GC 서열 | 15 | 66.67 | 2.50 |

예상값은 초기 코드의 정의에 따라 계산했다. 저GC 서열의 기록 거리 목록은 2, 1, 2이고 고GC 서열은 1, 4이다. 초기 코드의 출력 종명은 고정 문자열이며 합성 서열의 실제 종을 의미하지 않는다.

```bash
python - <<'PY'
import csv
with open('initial_summary.csv', newline='') as f:
    actual = list(csv.reader(f))
with open('expected_summary.csv', newline='') as f:
    expected = list(csv.reader(f))
assert actual == expected, (actual, expected)
print('PASS: synthetic example')
PY
```

이 검사는 초기 두 종 분석 코드의 작은 입력 예시를 확인하기 위한 것이며, 9종 확장 분석의 검증을 대신하지 않는다.
