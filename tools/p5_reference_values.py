# P5 (구조 동해석) 기준안의 해석해 유도와 참고값 계산.
# docs/design/03_verification_plan.md §1.5 의 "미리 계산한 참고값"을 이 스크립트로 다시 만들 수 있다.
# 판정에는 시험 코드(C++) 안에서 계산한 값을 쓴다. 이 스크립트는 유도 확인용이다.
#
# 필요: numpy, scipy, sympy (프로젝트 기본 환경에는 sympy, scipy 가 없을 수 있음)
#   python -m venv p5v && p5v\Scripts\python -m pip install numpy scipy sympy
#   p5v\Scripts\python tools\p5_reference_values.py
import numpy as np
import sympy as sp
from scipy import optimize, special


def axisym_stress(ur, uz, r, z, lam, mu):
    err, ezz, ett = sp.diff(ur, r), sp.diff(uz, z), ur / r
    grz = sp.diff(ur, z) + sp.diff(uz, r)
    tr = err + ezz + ett
    return (lam * tr + 2 * mu * err, lam * tr + 2 * mu * ett,
            lam * tr + 2 * mu * ezz, mu * grz)


def equilibrium(ur, uz, r, z, lam, mu):
    srr, stt, szz, srz = axisym_stress(ur, uz, r, z, lam, mu)
    f_r = sp.diff(srr, r) + sp.diff(srz, z) + (srr - stt) / r
    f_z = sp.diff(srz, r) + sp.diff(szz, z) + srz / r
    return sp.simplify(f_r), sp.simplify(f_z), (srr, stt, szz, srz)


def section(title):
    print("\n" + "=" * 70 + "\n" + title + "\n" + "=" * 70)


r, z, t = sp.symbols("r z t", real=True)
E, nu, rho, p = sp.symbols("E nu rho p", positive=True)
lam = E * nu / ((1 + nu) * (1 - 2 * nu))
mu = E / (2 * (1 + nu))

# ---------------------------------------------------------------- 4-A(a) 패치 시험
section("4-A(a) patch field: u_r = c1 r, u_z = c0 + c2 z  (body force 0?)")
c0, c1, c2 = sp.symbols("c0 c1 c2")
fr, fz, s = equilibrium(c1 * r, c0 + c2 * z, r, z, lam, mu)
print("residual body force:", fr, fz)
print("stresses (constant):", [sp.simplify(x) for x in s])

# ---------------------------------------------------------------- 4-B Lamé
section("4-B Lame thick cylinder, plane strain, inner pressure p")
a, b = sp.symbols("a b", positive=True)
ur = (1 + nu) * p * a**2 / (E * (b**2 - a**2)) * ((1 - 2 * nu) * r + b**2 / r)
fr, fz, (srr, stt, szz, srz) = equilibrium(ur, sp.Integer(0), r, z, lam, mu)
print("equilibrium residual:", fr, fz)
print("sigma_rr(a) =", sp.simplify(srr.subs(r, a)), "  sigma_rr(b) =", sp.simplify(srr.subs(r, b)))
print("sigma_tt =", sp.simplify(stt))
print("sigma_zz =", sp.simplify(szz))
vals = {E: 1, nu: sp.Rational(3, 10), p: 1, a: sp.Rational(1, 2), b: 1}
print("E=1, nu=0.3, p=1, a=0.5, b=1: u_r(a) = %.10f, sigma_tt(a) = %.10f" %
      (float(ur.subs(vals).subs(r, sp.Rational(1, 2))), float(stt.subs(vals).subs(r, sp.Rational(1, 2)))))

# ---------------------------------------------------------------- 4-J 균일 가속 막대
section("4-J(a) rod under uniform body force b_z = -rho*g (top z=L held, bottom z=0 free)")
g, L = sp.symbols("g L", positive=True)
ur = -nu * rho * g * z * r / E
uz = rho * g * z**2 / (2 * E) + nu * rho * g * r**2 / (2 * E)
srr, stt, szz, srz = axisym_stress(ur, uz, r, z, lam, mu)
f_r = sp.simplify(sp.diff(srr, r) + sp.diff(srz, z) + (srr - stt) / r)
f_z = sp.simplify(sp.diff(srz, r) + sp.diff(szz, z) + srz / r - rho * g)
print("residual (div sigma + b) with body force -rho g:", f_r, f_z)
print("stresses:", [sp.simplify(x) for x in (srr, stt, szz, srz)])

section("4-J(b) free rod, end pressure q on z=0 face, inertia relief")
q = sp.symbols("q", positive=True)
# 강체 가속도 A = q/(rho L) (+z), 상대 좌표계 체적력 -rho A. 응력 sigma_zz = -q (1 - z/L).
A_ = q / (rho * L)
szz_exact = -q * (1 - z / L)
ur = -nu * r * szz_exact / E
uz = sp.integrate(szz_exact / E, z) + nu * q * r**2 / (2 * E * L)
srr, stt, szz, srz = axisym_stress(ur, uz, r, z, lam, mu)
f_r = sp.simplify(sp.diff(srr, r) + sp.diff(srz, z) + (srr - stt) / r)
f_z = sp.simplify(sp.diff(srz, r) + sp.diff(szz, z) + srz / r - rho * A_)
print("residual with body force -rho A:", f_r, f_z)
print("stresses:", [sp.simplify(x) for x in (srr, stt, szz, srz)])
print("u_r =", sp.simplify(ur), "   u_z =", sp.simplify(uz))

# ---------------------------------------------------------------- 4-C 원판 굽힘 다항식 해
section("4-C circular plate, uniform pressure p on top, exact polynomial elasticity solution")
h = sp.symbols("h", positive=True)
deg = 7
coef = []
urp, uzp = 0, 0
for i in range(0, deg):
    for j in range(0, deg):
        if 2 * i + 1 + j <= deg:
            c = sp.Symbol(f"a_{i}_{j}")
            coef.append(c)
            urp += c * r ** (2 * i + 1) * z**j
        if 2 * i + j <= deg:
            c = sp.Symbol(f"b_{i}_{j}")
            coef.append(c)
            uzp += c * r ** (2 * i) * z**j
nuv = sp.Rational(3, 10)
lam_n, mu_n = lam.subs({E: 1, nu: nuv}), mu.subs({E: 1, nu: nuv})
srr, stt, szz, srz = axisym_stress(urp, uzp, r, z, lam_n, mu_n)
eqs = []
for expr in (sp.diff(srr, r) + sp.diff(srz, z) + (srr - stt) / r,
             sp.diff(srz, r) + sp.diff(szz, z) + srz / r):
    eqs += sp.Poly(sp.expand(expr * r), r, z).coeffs()
for zz, pz in ((h / 2, -1), (-h / 2, 0)):
    eqs += sp.Poly(sp.expand(szz.subs(z, zz) - pz), r).coeffs()
    eqs += sp.Poly(sp.expand(srz.subs(z, zz)), r).coeffs()
sol = sp.solve(eqs, coef, dict=True)[0]
urs, uzs = sp.expand(urp.subs(sol)), sp.expand(uzp.subs(sol))
free = sorted((urs.free_symbols | uzs.free_symbols) - {r, z, h}, key=str)
print("free constants after equilibrium + face tractions:", free)
# 남은 상수: 강체 이동, 반경 늘림, 순수 굽힘 등. r = a = 1 에서 중립면 고정(Kirchhoff 고정단 의미)으로 정한다.
fix = [uzs.subs({r: 1, z: 0}), sp.diff(uzs, r).subs({r: 1, z: 0}), urs.subs({r: 1, z: 0})]
sol2 = sp.solve(fix, free, dict=True)[0] if free else {}
urs, uzs = sp.expand(urs.subs(sol2)), sp.expand(uzs.subs(sol2))
rest = sorted((urs.free_symbols | uzs.free_symbols) - {r, z, h}, key=str)
urs, uzs = urs.subs({s_: 0 for s_ in rest}), uzs.subs({s_: 0 for s_ in rest})
print("remaining constants set to 0:", rest)
print("u_r(r,z) =", sp.collect(urs, [r, z]))
print("u_z(r,z) =", sp.collect(uzs, [r, z]))
# 확인
srr, stt, szz, srz = axisym_stress(urs, uzs, r, z, lam_n, mu_n)
chk = [sp.simplify(sp.diff(srr, r) + sp.diff(srz, z) + (srr - stt) / r),
       sp.simplify(sp.diff(srz, r) + sp.diff(szz, z) + srz / r),
       sp.simplify(szz.subs(z, h / 2) + 1), sp.simplify(szz.subs(z, -h / 2)),
       sp.simplify(srz.subs(z, h / 2)), sp.simplify(srz.subs(z, -h / 2))]
print("check (all 0):", chk)
for hv in (sp.Rational(1, 10), sp.Rational(1, 50)):
    D = hv**3 / (12 * (1 - nuv**2))
    w0 = -uzs.subs({r: 0, z: 0, h: hv})
    print("h/a = %s: w(0,0) = %.10e, Kirchhoff p a^4/(64D) = %.10e, ratio = %.6f, sigma_rr(0,-h/2) = %.6e" %
          (hv, float(w0), float(1 / (64 * D)), float(w0 * 64 * D), float(sp.limit(srr.subs({z: -hv / 2, h: hv}), r, 0))))

# ---------------------------------------------------------------- 4-E 고유진동수
section("4-E(a) rod axial modes, nu = 0, fixed z=0, free z=L: omega_n = (2n-1) pi c / (2L)")
for n in (1, 2, 3):
    print("n=%d: omega L / c = %.12f" % (n, (2 * n - 1) * np.pi / 2))

section("4-E(b) hollow cylinder plane-strain radial modes (inner a, outer b, traction free)")
def radial_det(k, a_, b_, nu_):
    lam_ = nu_ / ((1 + nu_) * (1 - 2 * nu_))
    mu_ = 1 / (2 * (1 + nu_))
    c11 = lam_ + 2 * mu_
    def srr(f, fp, rr):
        return c11 * k * fp(k * rr) + lam_ * f(k * rr) / rr
    J, Jp = special.j1, lambda x: special.jvp(1, x)
    Y, Yp = special.y1, lambda x: special.yvp(1, x)
    return srr(J, Jp, a_) * srr(Y, Yp, b_) - srr(Y, Yp, a_) * srr(J, Jp, b_)

for (a_, b_) in ((0.5, 1.0), (0.9, 1.0)):
    nu_ = 0.3
    c_p = np.sqrt((nu_ / ((1 + nu_) * (1 - 2 * nu_)) + 1 / (1 + nu_)))  # sqrt((lam+2mu)/rho), E=rho=1
    ks = np.linspace(0.05, 40, 40000)
    d = [radial_det(k, a_, b_, nu_) for k in ks]
    roots = []
    for i in range(len(ks) - 1):
        if np.sign(d[i]) != np.sign(d[i + 1]):
            roots.append(optimize.brentq(radial_det, ks[i], ks[i + 1], args=(a_, b_, nu_), xtol=1e-15))
        if len(roots) == 3:
            break
    print("a=%.2f b=%.2f nu=0.3 E=rho=1: omega = k*c_p =" % (a_, b_),
          ", ".join("%.10f" % (k * c_p) for k in roots))
    # 얇은 고리 근사 (확인용): omega^2 = E/(rho(1-nu^2) R^2)
    R = 0.5 * (a_ + b_)
    print("   thin-ring check omega ~ %.6f" % np.sqrt(1 / ((1 - nu_**2) * R**2)))

section("4-E(c) clamped Kirchhoff plate: J0(l) I1(l) + I0(l) J1(l) = 0  (report only)")
f = lambda x: special.j0(x) * special.i1(x) + special.i0(x) * special.j1(x)
for lo, hi in ((2.5, 4.0), (5.5, 7.0), (8.5, 10.0)):
    lr = optimize.brentq(f, lo, hi, xtol=1e-15)
    print("lambda = %.10f, lambda^2 = %.10f  (omega = lambda^2/a^2 sqrt(D/(rho h)))" % (lr, lr**2))

# ---------------------------------------------------------------- 4-H 시간 적분 이산 정확해
section("4-H Newmark average acceleration and HHT-alpha on a single mode (derived here)")


def hht_matrix(Om, al, xi=0.0):
    # 단자유도 u'' + 2 xi w u' + w^2 u = 0, w = 1, Om = w dt.
    # HHT-alpha (alpha in [-1/3, 0]):
    #   M a_{n+1} + (1+al) (C v_{n+1} + K u_{n+1}) - al (C v_n + K u_n) = 0
    #   u_{n+1} = u_n + dt v_n + dt^2 ((1/2 - beta) a_n + beta a_{n+1})
    #   v_{n+1} = v_n + dt ((1 - gamma) a_n + gamma a_{n+1})
    # gamma = (1 - 2 al)/2, beta = (1 - al)^2 / 4. al = 0 -> Newmark 평균가속도.
    dt = Om
    ga, be = (1 - 2 * al) / 2, (1 - al) ** 2 / 4
    K, C = 1.0, 2 * xi
    # 미지수 x = (u1, v1, a1), 상태 (u0, v0, a0)
    Aimp = np.array([[1, 0, -be * dt**2], [0, 1, -ga * dt], [(1 + al) * K, (1 + al) * C, 1]])
    Bexp = np.array([[1, dt, (0.5 - be) * dt**2], [0, 1, (1 - ga) * dt], [al * K, al * C, 0]])
    return np.linalg.solve(Aimp, Bexp)


for Om in (0.1, 1.0, 10.0):
    lamb = np.linalg.eigvals(hht_matrix(Om, 0.0))
    big = lamb[np.argsort(-abs(lamb))[:2]]
    wbar = abs(np.angle(big[0])) / Om
    print("Newmark Om=%5.2f: |lambda| = %.15f, wbar/w = %.12f, closed form 2 atan(Om/2)/Om = %.12f"
          % (Om, abs(big[0]), wbar, 2 * np.arctan(Om / 2) / Om))
for al in (-0.1, -0.3):
    rho_inf = max(abs(np.linalg.eigvals(hht_matrix(1e6, al))))
    print("HHT alpha=%.2f: spectral radius at Om=1e6 = %.8f, (1+al)/(1-al) = %.8f; at Om=0.1: %.12f"
          % (al, rho_inf, (1 + al) / (1 - al), max(abs(np.linalg.eigvals(hht_matrix(0.1, al))))))
