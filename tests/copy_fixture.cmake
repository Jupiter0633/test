# 시험용 프로젝트 루트를 새로 복사한다. 입력: SRC, DST
file(REMOVE_RECURSE "${DST}")
file(COPY "${SRC}/" DESTINATION "${DST}")
