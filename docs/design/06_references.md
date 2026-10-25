# 06. 참고문헌과 확인 수준

확인일: 2026-10-21. 확인 수단: 웹 검색 결과(서지 데이터베이스, 출판사 페이지, 저장소 기록, 인용 목록).
**이번 작업에서 본문 전체를 직접 읽어 확인한 문헌은 없습니다.** 따라서 수치(표, 그래프 값)를 통과 기준으로 쓰기 전에 반드시 원문을 확인해야 합니다.
(2026-10-22 추가: R-21, R-30, R-31 은 해당 쪽의 원문 이미지를 직접 읽어 확인했습니다. 책 전체를 읽은 것은 아닙니다.)

확인 수준:
- **원문 확인**: 해당 내용이 실린 쪽의 원문을 직접 읽어 확인.
- **서지 확인**: 저자·제목·학술지·권·쪽 등을 두 곳 이상 또는 출판사/공식 기록에서 확인. 내용은 미확인.
- **초록 확인**: 서지 + 초록(또는 공식 요약)의 내용까지 확인.
- **2차 출처**: 내용을 다른 문헌의 인용·요약으로만 확인.
- **미확인**: 서지나 내용을 확인하지 못함. 근거로 사용 금지.

| ID | 문헌 | 확인 수준 | 본 문서에서의 용도 | 비고 |
|---|---|---|---|---|
| R-01 | Th. von Kármán (1929), "The impact on seaplane floats during landing," NACA TN-321 | 서지 확인 (NASA NTRS, Caltech 저장소) | 배경 (운동량·부가질량 기반 충격 이론) | |
| R-02 | H. Wagner (1932), "Über Stoß- und Gleitvorgänge an der Oberfläche von Flüssigkeiten," ZAMM 12(4):193–215 (영문 번역: NACA Translation No. 1366) | 서지 확인 | 배경, 3-A 근사해 | |
| R-03 | R. Zhao, O. Faltinsen (1993), "Water entry of two-dimensional bodies," J. Fluid Mech. 246:593–612 | 초록 확인 (데드라이즈 4°–81° 상사해, 비선형 BEM) | 3-A 비교 대상 | 상사해 수치표는 원문 확인 필요 |
| R-04 | O. M. Faltinsen, *Hydrodynamics of High-Speed Marine Vehicles*, Cambridge Univ. Press (2005; 일부 목록은 2006), "Slamming, Whipping, and Springing" 장 | 서지 확인 | 배경 | 장 번호·출판연도 표기가 목록마다 다름 |
| R-05 | S.-L. Chuang (1966), "Experiments on flat-bottom slamming," J. Ship Res. 10 | 서지 일부 확인 (쪽수 10–17 / 10–27 로 출처마다 불일치), 내용은 2차 출처 | 공기 완충이 첨두압력을 낮추고 충격시간을 늘린다는 정성 근거 | 쪽수 원문 확인 필요 |
| R-06 | G. Lewison, W. M. Maclean (1968), "On the cushioning of water impact by entrapped air," J. Ship Res. | 서지 일부 확인 (권·쪽 미확인), 요지 확인 (2D 조건 첨두압력이 기존 이론보다 훨씬 작음) | 공기 완충 정성 근거 | 권·쪽 미확인 |
| R-07 | C. W. Hirt, B. D. Nichols (1981), "Volume of fluid (VOF) method for the dynamics of free boundaries," J. Comput. Phys. 39(1):201–225, doi:10.1016/0021-9991(81)90145-5 | 초록 확인 | VOF 기반 | |
| R-08 | J. Roenby, H. Bredmose, H. Jasak (2016), "A computational method for sharp interface advection," R. Soc. Open Sci. 3:160405, doi:10.1098/rsos.160405 | 초록 확인 | 기하 VOF 후보 (isoAdvector) | |
| R-09 | R. Mittal, G. Iaccarino (2005), "Immersed boundary methods," Annu. Rev. Fluid Mech. 37:239–261, doi:10.1146/annurev.fluid.37.061903.175743 | 서지 확인 | IBM 분류 (연속 강제 vs 직접 강제/날카로운 경계) | 분류 내용은 2차 출처 |
| R-10 | C. S. Peskin (2002), "The immersed boundary method," Acta Numerica 11:479–517 | 초록 확인 | IBM 기원 (확산형) | |
| R-11 | Z. H. Ma, D. M. Causon, L. Qian, C. G. Mingham, T. Mai, D. Greaves, A. Raby (2016), "Pure and aerated water entry of a flat plate," Phys. Fluids 28(1), doi:10.1063/1.4940043 | 서지 확인, 내용 미확인 | 평판 입수 수치 해석 참조 후보 (3D 확장) | 공기 압축성 처리 내용 원문 확인 필요 |
| R-12 | S. Okada, Y. Sumi (2000), "On the water impact and elastic response of a flat plate at small impact angles," J. Mar. Sci. Technol. 5:31–39, doi:10.1007/s007730070019 | 초록 확인 (0°–4° 경사, 압력 이력 3유형: Wagner형/공기포획형/중간형, 최대 변형률은 각도에 덜 민감) | 경사 민감도, 5-E 정성 비교 | |
| R-13 | A. Korobkin (1996), "Acoustic approximation in the slamming problem," J. Fluid Mech. 318:165–188 | 초록 확인 | 물 압축성 초기 단계 한계 근거 | |
| R-14 | T. Mai, C. Mai, A. Raby, D. Greaves (2019), "Aeration effects on water-structure impacts: Part 1. Drop plate impacts," Ocean Eng. 193:106600 | 초록 확인 (정사각 평판, 데드라이즈 0, 낙하속도 약 1–7 m/s, 포기율 0–1.6%; 포기가 첨두하중을 크게 낮추고, 충격량 변동은 첨두하중보다 작음) | 충격량 지표 채택 근거, 7-B | 저자 표기가 저장소마다 다름 |
| R-15 | O. M. Faltinsen (1997), "The effect of hydroelasticity on ship slamming," Phil. Trans. R. Soc. A 355(1724):575–591, doi:10.1098/rsta.1997.0026 | 서지 확인, 내용은 학회 버전 요약(2차 출처) | 수탄성 시 응력이 입수속도에 선형 비례 | |
| R-16 | P. Causin, J.-F. Gerbeau, F. Nobile (2005), "Added-mass effect in the design of partitioned algorithms for fluid–structure problems," Comput. Methods Appl. Mech. Eng. 194(42–44):4506–4527, doi:10.1016/j.cma.2004.12.005 | 초록 확인 | 분할형 연성의 부가질량 불안정 (6-F, D-13) | |
| R-17 | J. C. Martin, W. J. Moyce (1952), "Part IV. An experimental study of the collapse of liquid columns on a rigid horizontal plane," Phil. Trans. R. Soc. A 244:312–324 | 초록 확인 (2D 직사각·반원 단면 붕괴, **원기둥의 축대칭 붕괴** 포함) | 2-D, 2-E | 측정 불확도 제시 여부 미확인 |
| R-18 | U. Ghia, K. N. Ghia, C. T. Shin (1982), "High-Re solutions for incompressible flow using the Navier–Stokes equations and a multigrid method," J. Comput. Phys. 48(3):387–411 | 초록 확인 | 1-C | 표 수치는 원문 확인 필요 |
| R-19 | A. Prosperetti (1981), "Motion of two superposed viscous fluids," Phys. Fluids 24(7):1217–1223, doi:10.1063/1.863522 | 서지 확인, 내용은 2차 출처 (라플라스 공간 해) | 2-C 감쇠 정성 비교 | 평면 문제. 축대칭 적용 미확인 |
| R-20 | A. W. Leissa (1969), *Vibration of Plates*, NASA SP-160 | 서지 확인 (NASA NTRS) | 4-B | 원판 진동수 표 원문 확인 필요 |
| R-21 | H. Lamb, *Hydrodynamics*, 6th ed., Cambridge Univ. Press (1932). 부가질량: 구 Art. 92 p.124 식 (3)–(4) (밀어낸 질량의 절반), 얇은 원판 Art. 102 p.139 식 (20) (2T = (8/3)ρa³U²), 타원체 Art. 114 p.153 식 (6)–(12) (k = α₀/(2−α₀), 구에서 1/2) | **원문 확인** (2026-10-22, 인터넷 아카이브 스캔본 `dli.ernet.3284` 의 해당 쪽 이미지를 직접 읽음. 1932년판은 공개 자료). 이전 수준: 2차 출처 | 1-E, 1-F | 납작한 회전타원체가 원판 값으로 수렴함을 수치 적분으로 확인 (두께비 0.01 에서 99.7 %, 제 계산). 자유수면 위 원판의 "절반" 값은 **여전히 미확인** |
| R-22 | J. W. Glasheen, T. A. McMahon (1996), "Vertical water entry of disks at low Froude numbers," Phys. Fluids 8(8):2078–2083 | 서지 확인, 내용 미확인 | 3-H, 6-G 후보 | |
| R-23 | E. V. Ermanyuk, M. Ohkusu (2005), "Impact of a disk on shallow water," J. Fluids Struct. 20(3):345–357 | 서지 확인, 내용 미확인 | 3-H, 6-G 후보 (얕은 수심) | |
| R-24 | H. C. Mayer, R. Krechetnikov (2018), "Flat plate impact on water," J. Fluid Mech. 850 | 서지 확인 (권 확인, 쪽 미확인), 초록 일부만 확인 | 3-B, 3-H 후보 | 판 형상(원판 여부) 미확인 |
| R-25 | I. B. Celik et al. (2008), "Procedure for estimation and reporting of uncertainty due to discretization in CFD applications," J. Fluids Eng. 130(7):078001, doi:10.1115/1.2960953 | 서지 확인 (2차 출처) | GCI 절차 | |
| R-26 | P. J. Roache (1994), "Perspective: A method for uniform reporting of grid refinement studies," J. Fluids Eng. | 서지 일부 확인 (권·쪽 미확인) | GCI 원전 | |
| R-27 | P. J. Roache (2002), "Code verification by the method of manufactured solutions," J. Fluids Eng. 124 | 서지 일부 확인 (포럼 BibTeX — 쪽 미확인) | MMS (1-A) | |
| R-28 | ASME V&V 20-2009 (R2016), *Standard for Verification and Validation in Computational Fluid Dynamics and Heat Transfer* | 공식 요약 확인 (ASME 보도자료, ANSI) | 검증 불확도 u_val 기반 판정 틀 | 표준 원문 미확인. 판정식 세부(포함계수 등)는 원문 확인 후 확정 |
| R-30 | H. Lamb (1932), R-21 과 같은 책. Stokes 항력과 유속장: Art. 337 p.598 식 (10)·(11)·(15) (6πμaU, 바깥 벽이 지름의 여러 배 떨어져도 결과가 크게 달라진다는 주의). Oseen 보정 6πμUa(1 + 3/8·Ua/ν): Art. 343 p.617 식 (16) | **원문 확인** (2026-10-22, R-21 과 같은 스캔본) | 1-G | |
| R-31 | H. Lamb (1932), R-21 과 같은 책. 원형 수조 대칭 모드 J₀′(ka) = 0, ka = 3.832: Art. 191 p.285 식 (6)–(7a). 유한 깊이 σ² = gk tanh kh, 같은 k: Art. 257 p.440–441 식 (4). 위아래가 막힌 두 유체 σ² = gk(ρ−ρ′)/(ρ coth kh + ρ′ coth kh′): Art. 231 p.371 식 (11) | **원문 확인** (2026-10-22, R-21 과 같은 스캔본) | 2-C | 두 식의 조합은 유도 (`03` §1.3) |
| R-32 | G. D. Weymouth, D. K.-P. Yue (2010), "Conservative Volume-of-Fluid method for free-surface simulations on Cartesian-grids," J. Comput. Phys. 229(8):2853–2865 | **미확인** (서지는 기억에 의존, 원문·서지 데이터베이스 미확인) | P3 VOF 수송 방식 (방향 분할 + c·∂u/∂x 보정)의 착상 | 기준 수치로 쓰지 않음. 보존·범위는 `03` §1.3 시험(2-B)과 단위 시험으로 직접 확인 |
| R-33 | B. van Leer (1974), "Towards the ultimate conservative difference scheme. II. Monotonicity and conservation combined in a second-order scheme," J. Comput. Phys. 14(4):361–370 | **미확인** (서지는 기억에 의존) | P3 제한자(limiter) 대류 방식 (2-H, 2-D) 의 착상 | 기준 수치로 쓰지 않음. 정확도는 2-H 로 직접 확인 |
| R-34 | M. Rudman (1998), "A volume-tracking method for incompressible multifluid flows with large density variations," Int. J. Numer. Methods Fluids 28(2):357–378 | **미확인** (서지는 기억에 의존) | P3 일관 질량–운동량 수송 (방식 B, 필요할 때만) 의 착상 | 기준 수치로 쓰지 않음 |
| R-35 | U. Jain, P. Vega-Martínez, D. van der Meer (2021), "Air entrapment and its effect on pressure impulses in the slamming of a flat disc on water," J. Fluid Mech. 928, A31, doi:10.1017/jfm.2021.846 (CC BY) | 본문 일부 확인 (출판사 공개 HTML 의 인용 문장, 요약 도구 경유, 2026-10-23). 그림·표 수치 미확인 | 원판 입수의 공기막 두께 (h₀ ≈ 7×10⁻³ R), 0.5 m/s 이상에서 공기막이 뚫리거나 밀려 나감, 중심 첨두압 ≈ 100 ρV², 힘 충격량 ≈ 고전 부가질량 (4/3)ρR³V | 3-H 실험 비교 후보 (반지름 4 cm, 0–3 m/s, 우리 조건과 가까움). 정량 비교 전 원문 PDF 의 수치·불확도 확인 필요 |
| R-36 | Shi 외 (2025), "Mechanism of the sharp rise and fall of the impact force during water impact of a flat plate: compression-expansion behavior of trapped air layer," J. Fluids Struct. 138, 104393 | 서지·초록만 | 갇힌 공기층의 압축·팽창이 초기 힘의 급상승·급하강을 만듦 (수치) | 정성 근거만 |
| R-37 | Kim, Kim, Kim (2021), J. Fluid Mech. — 표면 무늬가 있는 원판의 입수 | 2차 출처 (검색 요약) | 갇힌 공기 부피가 클수록 첨두 낮고 주기 길어짐 | 서지 확인 필요 |
| R-38 | H. Sun, W. Ding, X. Zhao, Z. Sun (2022), "Numerical Study of Flat Plate Impact on Water Using a Compressible CIP–IBM–Based Model," J. Mar. Sci. Eng. 10(10), 1462, doi:10.3390/jmse10101462 (CC BY) | 본문 확인 (출판사 HTML, 2026-10-26) | 평판 입수 수치 해석의 격자 (2 mm, 반폭/Δ = 100) 와 Δt (5×10⁻⁶ s), 공기 쿠션과 압력 진동의 관계 | 격자·Δt 설정의 문헌 사례 (기준 수치로 쓰지 않음) |
| R-29 | Shi, Qu, Liu, Zheng, Li (2026), "Influence of plate size on air cushioning during water entry," Phys. Rev. Fluids 11:104801 | 검색 결과로만 존재 확인, 내용 미확인 | 7-E 후보 | 최근 문헌. 원문 확인 전 사용 금지 |

## 확인하지 못해 본문에서 근거로 쓰지 않은 내용

- 축대칭 기하 VOF 의 r-가중 재구성 기법 문헌
- 고밀도비 일관 질량–운동량 수송 기법 문헌
- 원통 쉘 축대칭 고유진동수 해석해 출처 (4-C)
- 습윤 고유진동수 평가 출처 (5-F)
- 규칙파 이론(Stokes 2차 등) 출처 (7-C)
- 단방향 연성 유효성 임계값의 문헌 근거 (5-D)
- Verhagen (1967) 평판 충격 논문 — 검색으로 서지 확인 실패
