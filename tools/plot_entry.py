"""개발용 그림: hcw_dev entry-snap 이 쓴 snap_NNN.csv, history.csv 로 입수 과정 그림을 만든다.

사용법:  python tools/plot_entry.py <snapshot 폴더> [r_max/a] [z_min/a] [z_max/a] [t_from ms] [t_to ms]
  t_from, t_to 를 주면 그 구간의 스냅샷만 그린다 (최대 16 장).
출력 (같은 폴더):
  water.png    시각별 물 (파랑) – 공기 (흰색) 분포, 물체 (회색), 속도 화살표
  pressure.png 시각별 게이지 압력 색 지도 (빨강 +, 파랑 −), 물–공기 경계선
  history.png  물체 힘과 바닥 평균·포켓 압력의 시간 이력 (그림 시각 표시)
  diag.png     최저 절대압력, 압력 하한 발동 누적, 축 위 공기층 두께, 포켓 부피
그림 안 글자는 영어 (프로그램 출력 규칙).
"""

import glob
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


def read_snapshot(path):
    with open(path, encoding="utf-8") as f:
        head = f.readline().lstrip("# ").split()
    meta = {k: float(v) for k, v in (item.split("=") for item in head)}
    data = np.loadtxt(path, delimiter=",", skiprows=2)
    nr, nz = int(meta["nr"]), int(meta["nz"])
    fields = {name: data[:, k].reshape(nr, nz) for k, name in
              enumerate(["r", "z", "phi", "water", "p", "u", "w", "active"])}
    return meta, fields


def body_patch(meta):
    a, zc, hh = meta["a"], meta["zc"], meta["half_height"]
    return Rectangle((0.0, (zc - hh) / a), 1.0, 2.0 * hh / a, facecolor="0.55", edgecolor="k", lw=0.8, zorder=5)


def main():
    folder = sys.argv[1]
    r_max = float(sys.argv[2]) if len(sys.argv) > 2 else 1.6
    z_min = float(sys.argv[3]) if len(sys.argv) > 3 else -0.6
    z_max = float(sys.argv[4]) if len(sys.argv) > 4 else 0.6
    t_from = float(sys.argv[5]) if len(sys.argv) > 5 else -1.0
    t_to = float(sys.argv[6]) if len(sys.argv) > 6 else 1e9
    snaps = sorted(glob.glob(os.path.join(folder, "snap_*.csv")))
    if not snaps:
        sys.exit("no snap_*.csv in " + folder)
    loaded = [read_snapshot(p) for p in snaps]
    loaded = [x for x in loaded if t_from - 1e-9 <= x[0]["t"] * 1e3 <= t_to + 1e-9]
    if len(loaded) > 16:
        stride = math.ceil(len(loaded) / 16)
        loaded = loaded[::stride]
    a = loaded[0][0]["a"]
    cols = 4
    rows = math.ceil(len(loaded) / cols)

    # 압력 색 범위: 보이는 영역 안의 최대 크기 (모든 시각 공통)
    pmax = 1.0
    for meta, f in loaded:
        sel = (f["r"] / a <= r_max) & (f["z"] / a >= z_min) & (f["z"] / a <= z_max) & (f["active"] > 0)
        if sel.any():
            pmax = max(pmax, np.abs(f["p"][sel]).max())

    for kind in ("water", "pressure"):
        fig, axes = plt.subplots(rows, cols, figsize=(4.0 * cols, 3.0 * rows), squeeze=False)
        for ax in axes.flat:
            ax.set_visible(False)
        for ax, (meta, f) in zip(axes.flat, loaded):
            ax.set_visible(True)
            r, z = f["r"][:, 0] / a, f["z"][0, :] / a
            # 늘어나는 격자도 그리도록 셀 중심 좌표로 칸을 칠한다 (shading="nearest")
            water = np.where(f["phi"] > 0.5, f["water"], np.nan)
            if kind == "water":
                ax.pcolormesh(r, z, water.T, cmap="Blues", vmin=0, vmax=1, shading="nearest")
                ax.set_aspect("equal")
                # 속도 화살표: 칸 크기와 상관없이 0.12a 간격에 가장 가까운 셀
                ii = np.unique(np.clip(np.searchsorted(r, np.arange(r[0], r_max, 0.12)), 0, len(r) - 1))
                jj = np.unique(np.clip(np.searchsorted(z, np.arange(z_min, z_max, 0.12)), 0, len(z) - 1))
                R, Z = np.meshgrid(r[ii], z[jj], indexing="ij")
                U, W = f["u"][np.ix_(ii, jj)], f["w"][np.ix_(ii, jj)]
                keep = (R <= r_max) & (Z >= z_min) & (Z <= z_max) & (f["phi"][np.ix_(ii, jj)] > 0.5)
                ax.quiver(R[keep], Z[keep], U[keep], W[keep], color="tab:red", scale=60, width=0.004, zorder=6)
            else:
                p = np.where((f["active"] > 0) & (f["phi"] > 0.5), f["p"], np.nan) / 1e3
                im = ax.pcolormesh(r, z, p.T, cmap="RdBu_r", shading="nearest",
                                   norm=TwoSlopeNorm(vcenter=0.0, vmin=-pmax / 1e3, vmax=pmax / 1e3))
                ax.set_aspect("equal")
                ax.contour(r, z, water.T, levels=[0.5], colors="k", linewidths=0.8)
            ax.add_patch(body_patch(meta))
            ax.set_xlim(0, r_max)
            ax.set_ylim(z_min, z_max)
            ax.set_title(f"t = {meta['t'] * 1e3:.2f} ms", fontsize=10)
            ax.tick_params(labelsize=7)
        for ax in axes[-1]:
            ax.set_xlabel("r / a", fontsize=8)
        for ax in axes[:, 0]:
            ax.set_ylabel("z / a", fontsize=8)
        if kind == "pressure":
            fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.6, label="gauge pressure (kPa)")
            fig.suptitle("Gauge pressure (black line: water-air interface, grey: disk)", fontsize=12)
        else:
            fig.suptitle("Water fraction (blue = water, white = air, grey = disk, red arrows = velocity)",
                         fontsize=12)
        out = os.path.join(folder, kind + ".png")
        fig.savefig(out, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out)

    hist = np.loadtxt(os.path.join(folder, "history.csv"), delimiter=",", skiprows=1)
    fig, ax1 = plt.subplots(figsize=(9, 4))
    ax1.plot(hist[:, 0] * 1e3, hist[:, 1], "k-", lw=1.2, label="force on disk (N)")
    ax1.set_xlabel("t (ms)")
    ax1.set_ylabel("force (N)")
    ax2 = ax1.twinx()
    ax2.plot(hist[:, 0] * 1e3, hist[:, 2] / 1e3, "tab:red", lw=1.0, label="bottom mean pressure (kPa)")
    ax2.plot(hist[:, 0] * 1e3, hist[:, 3] / 1e3, "tab:blue", lw=1.0, ls="--", label="trapped pocket pressure (kPa)")
    ax2.set_ylabel("gauge pressure (kPa)")
    for meta, _ in loaded:
        ax1.axvline(meta["t"] * 1e3, color="0.85", lw=0.6, zorder=0)
    lines = ax1.get_lines()[:1] + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], loc="upper right", fontsize=8)
    ax1.set_title("Disk entry history (grey lines: snapshot times)")
    out = os.path.join(folder, "history.png")
    fig.savefig(out, dpi=110, bbox_inches="tight")
    plt.close(fig)
    print("wrote", out)

    if hist.shape[1] >= 8:
        t = hist[:, 0] * 1e3
        fig, axs = plt.subplots(4, 1, figsize=(9, 9), sharex=True)
        axs[0].plot(t, hist[:, 1], "k-", lw=1)
        axs[0].set_ylabel("force (N)")
        axs[1].plot(t, hist[:, 5] / 1e3, "tab:purple", lw=1)
        axs[1].axhline(1.0, color="tab:red", ls=":", lw=1)
        axs[1].set_ylabel("min abs. pressure (kPa)")
        axs[1].text(0.01, 0.05, "red dotted: 1 kPa air pressure floor", transform=axs[1].transAxes, fontsize=8)
        axs[2].plot(t, np.diff(np.concatenate([[0.0], hist[:, 6]])), "tab:orange", lw=1)
        axs[2].set_ylabel("floored cells per step")
        axs[3].plot(t, hist[:, 7] / a, "tab:blue", lw=1, label="air gap on axis / a")
        axs[3].plot(t, hist[:, 4] * 2 * np.pi / a**3, "tab:green", lw=1, label="pocket volume / a^3")
        axs[3].set_ylabel("gap, volume")
        axs[3].legend(fontsize=8)
        axs[3].set_xlabel("t (ms)")
        for ax in axs:
            ax.grid(alpha=0.3)
        fig.suptitle("Entry diagnostics")
        out = os.path.join(folder, "diag.png")
        fig.savefig(out, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print("wrote", out)


if __name__ == "__main__":
    main()
