# 분석 결과

최종 9종 유전체의 분석 결과를 모은 폴더이다. 새 분석을 수행할 때는 별도의 출력 폴더를 지정한다.

| 파일 | 내용 |
|---|---|
| [strain_summary.csv](strain_summary.csv) | 균주별 GC 함량, 평균 종결 간격, 셔플 비교, 코딩 밀도 |
| [frames.csv](frames.csv) | 관측 및 셔플 서열의 6-frame 지표 |
| [shuffle.csv](shuffle.csv) | 셔플 반복별 통합 지표 |
| [regions.csv](regions.csv) | CDS 및 비CDS 구간 지표 |
| [shifts.csv](shifts.csv) | 검증 CDS의 읽기 틀 이동 분석 |
| [genes.csv](genes.csv) | CDS별 길이 및 검증 상태 |
| [correlations.csv](correlations.csv) | Spearman 순위상관계수 |
| [lifestyle.csv](lifestyle.csv) | 생활 방식 그룹별 요약 |
| [run_info.json](run_info.json) | 입력 accession·체크섬 및 실행 조건 |

## 그래프 생성



```bash
python plot.py results
```

결과 그래프 8개는 `results/figures/`에 생성된다. 
