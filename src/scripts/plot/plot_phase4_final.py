"""
Generate Phase 4 adaptive kernel charts for experiment_report.md.
Data from remote scripts/run/sweep_threads.ps1 on 2026-06-05, all MD5 PASS.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


threads = [1, 2, 4, 8, 16, 32, 64]

# test.fasta, threshold=0.85, one remote run per thread count.
p4_times_t1 = [0.080, 0.060, 0.040, 0.040, 0.030, 0.040, 0.040]

# test2.fasta, average across thresholds 0.80/0.85/0.90/0.95.
p4_times_t2 = [65.962, 34.812, 21.228, 12.845, 8.295, 5.620, 5.778]


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
OUTDIR = os.path.join(REPO_ROOT, "report_images")
os.makedirs(OUTDIR, exist_ok=True)


def speedup(times):
    return [times[0] / value for value in times]


def efficiency(times):
    sp = speedup(times)
    return [sp[i] / threads[i] * 100.0 for i in range(len(threads))]


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
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


def plot_speedup(times, title, filename):
    sp = speedup(times)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(threads, sp, "o-", color=colors["speedup"], label="阶段四实测加速比", linewidth=2.2, markersize=7)
    for th, value in zip(threads, sp):
        ax.annotate(f"{value:.2f}x", (th, value), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8)
    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("加速比 (Speedup)")
    ax.set_title(title)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(threads)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25, which="both")
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


def plot_efficiency(times, title, filename):
    eff = efficiency(times)
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(x, eff, 0.55, color=colors["efficiency"], edgecolor="white", linewidth=0.6)
    ax.axhline(y=100, color=colors["ideal"], linestyle="--", linewidth=1.2, label="理想效率 (100%)")
    for bar in bars:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + 2, f"{h:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("并行效率 (%)")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, 120)
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


plot_time(p4_times_t1, "阶段四执行时间 vs 线程数 (test.fasta, threshold=0.85)", "p4_exec_time_t1.png", "{:.3f}s")
plot_speedup(p4_times_t1, "阶段四加速比曲线 (test.fasta)", "p4_speedup_t1.png")
plot_efficiency(p4_times_t1, "阶段四并行效率 (test.fasta)", "p4_efficiency_t1.png")

plot_time(p4_times_t2, "阶段四执行时间 vs 线程数 (test2.fasta, 四阈值平均)", "p4_exec_time_t2.png", "{:.2f}s")
plot_speedup(p4_times_t2, "阶段四加速比曲线 (test2.fasta, 四阈值平均)", "p4_speedup_t2.png")
plot_efficiency(p4_times_t2, "阶段四并行效率 (test2.fasta)", "p4_efficiency_t2.png")

print("All Phase 4 charts saved successfully!")
