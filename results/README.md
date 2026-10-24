# 계산 결과 폴더

`hcw prepare` (이후 단계에서는 `hcw run`) 가 아래 구조로 결과 폴더를 만듭니다. 결과 파일은 용량 때문에 git 에 넣지 않습니다.

```
results/<dev|verification|paper>/<입력 이름>/<UTC 시각>/
  run_meta.txt   실행 상태, 입력 해시, 코드 커밋, 컴파일러, 기기 정보
  input.case     실행에 쓴 입력 파일 사본
```
