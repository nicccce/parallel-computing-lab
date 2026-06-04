"""
Generate optional Phase 3 postings performance charts for experiment_report.md.
Data from remote scripts/run/run_all.ps1 on 2026-06-04/05, all MD5 PASS.

Produces 8 charts:
  1. test.fasta  execution time
  2. test.fasta  speedup
  3. test.fasta  efficiency
  4. test2.fasta execution time
  5. test2.fasta speedup
  6. test2.fasta efficiency
  7. Phase 2 vs Phase 3 execution time comparison (test2.fasta)
  8. Phase 2 vs Phase 3 speedup comparison (test2.fasta)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Optional Phase 3 postings measured data from run_all.ps1.
threads = [1, 2, 4, 8, 16, 32, 64]

# test.fasta, threshold=0.85, one remote run per thread count.
p3_times_t1 = [0.080, 0.060, 0.050, 0.040, 0.040, 0.040, 0.040]

# test2.fasta, average across thresholds 0.80/0.85/0.90/0.95.
p3_times_t2 = [65.190, 34.565, 21.228, 12.712, 8.558, 5.502, 5.722]

# Phase 2 multi-threshold run_all.ps1 reference for apples-to-apples comparison.
# Data points are averages across thresholds 0.80/0.85/0.90/0.95, all MD5 PASS.
p2_times_t2 = [64.775, 34.820, 20.348, 13.000, 8.215, 5.338, 5.740]


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "report_images")

outdir = sys.argv[1] if len(sys.argv) >= 2 else DEFAULT_OUTDIR
os.makedirs(outdir, exist_ok=True)


def speedup(times):
    return [times[0] / t for t in times]


def efficiency(times):
    sp = speedup(times)
    return [sp[i] / threads[i] * 100 for i in range(len(threads))]


p3_sp_t1 = speedup(p3_times_t1)
p3_eff_t1 = efficiency(p3_times_t1)
p3_sp_t2 = speedup(p3_times_t2)
p3_eff_t2 = efficiency(p3_times_t2)
p2_sp_t2 = speedup(p2_times_t2)


plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans", "Arial"],
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

colors = {
    "time": "#2d7dd2",
    "speedup": "#27ae60",
    "ideal": "#95a5a6",
    "efficiency": "#f39c12",
    "phase2": "#bdc3c7",
    "phase3": "#e74c3c",
}

x = np.arange(len(threads))


def plot_time(times, title, filename, label_fmt):
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(x, times, 0.55, color=colors["time"], edgecolor="white", linewidth=0.6)
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + max(times) * 0.025,
            label_fmt.format(h),
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("执行时间 (s)")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(threads)
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, max(times) * 1.2)
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)


def plot_speedup(sp, title, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(
        threads,
        sp,
        "o-",
        color=colors["speedup"],
        label="可选 Phase 3 postings 实测加速比",
        linewidth=2.2,
        markersize=7,
    )
    for th, value in zip(threads, sp):
        ax.annotate(
            f"{value:.2f}×",
            (th, value),
            textcoords="offset points",
            xytext=(0, 8),
            ha="center",
            fontsize=8,
            color=colors["speedup"],
        )
    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("加速比 (Speedup)")
    ax.set_title(title)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(threads)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25, which="both")
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)


def plot_efficiency(eff, title, filename):
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(x, eff, 0.55, color=colors["efficiency"], edgecolor="white", linewidth=0.6)
    ax.axhline(y=100, color=colors["ideal"], linestyle="--", linewidth=1.2, label="理想效率 (100%)")
    for bar in bars:
        h = bar.get_height()
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            h + 2,
            f"{h:.1f}%",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("并行效率 (%)")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, 120)
    fig.savefig(os.path.join(outdir, filename))
    plt.close(fig)


plot_time(
    p3_times_t1,
    "可选 Phase 3 postings 执行时间 vs 线程数 (test.fasta, threshold=0.85)",
    "p3_exec_time_t1.png",
    "{:.3f}s",
)
plot_speedup(p3_sp_t1, "可选 Phase 3 postings 加速比曲线 (test.fasta)", "p3_speedup_t1.png")
plot_efficiency(p3_eff_t1, "可选 Phase 3 postings 并行效率 (test.fasta)", "p3_efficiency_t1.png")

plot_time(
    p3_times_t2,
    "可选 Phase 3 postings 执行时间 vs 线程数 (test2.fasta, 四阈值平均)",
    "p3_exec_time_t2.png",
    "{:.2f}s",
)
plot_speedup(p3_sp_t2, "可选 Phase 3 postings 加速比曲线 (test2.fasta, 四阈值平均)", "p3_speedup_t2.png")
plot_efficiency(p3_eff_t2, "可选 Phase 3 postings 并行效率 (test2.fasta, 四阈值平均)", "p3_efficiency_t2.png")

# Phase 2 vs Phase 3 execution time comparison on test2.fasta.
fig, ax = plt.subplots(figsize=(10, 6))
w = 0.35
bars2 = ax.bar(
    x - w / 2,
    p2_times_t2,
    w,
    color=colors["phase2"],
    edgecolor="white",
    linewidth=0.5,
    label="Phase 2",
)
bars3 = ax.bar(
    x + w / 2,
    p3_times_t2,
    w,
    color=colors["phase3"],
    edgecolor="white",
    linewidth=0.5,
    label="可选 Phase 3 postings",
)

for bar in bars3:
    h = bar.get_height()
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        h + max(p3_times_t2) * 0.018,
        f"{h:.2f}s",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#c0392b",
        weight="bold",
    )

ax.set_xlabel("线程数 (Thread Count)", fontsize=11, labelpad=10)
ax.set_ylabel("执行时间 (s)", fontsize=11, labelpad=10)
ax.set_title("Phase 2 vs 可选 Phase 3 postings 执行时间对比 (test2.fasta, 四阈值平均)", fontsize=14, pad=15, weight="bold")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.3)

axins = ax.inset_axes([0.45, 0.45, 0.4, 0.4])
axins.bar(x[-3:] - w / 2, p2_times_t2[-3:], w, color=colors["phase2"])
axins.bar(x[-3:] + w / 2, p3_times_t2[-3:], w, color=colors["phase3"])
axins.set_xticks(x[-3:])
axins.set_xticklabels(threads[-3:])
axins.set_title("16-64 线程放大", fontsize=10)
axins.grid(axis="y", alpha=0.3)
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.savefig(os.path.join(outdir, "p3_vs_p2_exec_time.png"))
plt.close(fig)

# Phase 2 vs Phase 3 speedup comparison on test2.fasta.
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.5)
ax.plot(threads, p2_sp_t2, "o--", color=colors["phase2"], label="Phase 2 加速比", linewidth=2)
ax.plot(
    threads,
    p3_sp_t2,
    "s-",
    color=colors["phase3"],
    label="可选 Phase 3 postings 加速比",
    linewidth=2.5,
    markersize=8,
)

for th, value in zip(threads, p3_sp_t2):
    ax.annotate(
        f"{value:.1f}×",
        (th, value),
        textcoords="offset points",
        xytext=(0, -15),
        ha="center",
        fontsize=9,
        weight="bold",
        color=colors["phase3"],
    )

ax.set_xlabel("线程数 (Thread Count)", fontsize=11, labelpad=10)
ax.set_ylabel("加速比 (Speedup)", fontsize=11, labelpad=10)
ax.set_title("Phase 2 vs 可选 Phase 3 postings 加速比对比 (test2.fasta, 四阈值平均)", fontsize=14, pad=15, weight="bold")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc="upper left")
ax.grid(True, alpha=0.3, which="both")
fig.savefig(os.path.join(outdir, "p3_vs_p2_speedup.png"))
plt.close(fig)

print("All Phase 3 charts saved successfully!")
