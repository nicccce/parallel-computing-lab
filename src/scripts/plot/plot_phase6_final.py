"""
Generate Phase 6 comprehensive tuning charts for experiment_report.md.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


threads = [1, 2, 4, 8, 16, 32, 64]

# Phase 6 data
p6_times_t1 = [0.05, 0.03, 0.02, 0.01, 0.01, 0.013, 0.013]
p6_times_t2 = [11.088, 5.785, 3.425, 1.862, 1.112, 0.810, 0.832]

# Phase 5 data for comparison
p5_times_t1 = [0.080, 0.050, 0.030, 0.020, 0.020, 0.020, 0.020]
p5_times_t2 = [63.493, 33.045, 20.008, 11.325, 7.042, 4.887, 8.380]

# Phase 1 baseline data
p1_times_t2 = [184.910, 105.410, 56.573, 30.967, 18.547, 11.587, 9.713]


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
    "time": "#9b59b6",
    "speedup": "#27ae60",
    "ideal": "#bdc3c7",
    "efficiency": "#f1c40f",
    "p5_time": "#2d7dd2",
    "p6_time": "#9b59b6",
    "p5_speedup": "#2d7dd2",
    "p6_speedup": "#9b59b6",
    "p1_time": "#e74c3c",
    "p1_speedup": "#e74c3c",
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
    ax.plot(threads, sp, "o-", color=colors["speedup"], label="阶段六实测加速比", linewidth=2.2, markersize=7)
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


def plot_comparison_time(times_old, times_new, label_old, label_new, title, filename):
    fig, ax = plt.subplots(figsize=(10, 5))
    w = 0.35
    bars1 = ax.bar(x - w/2, times_old, w, label=label_old, color=colors["p5_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)
    bars2 = ax.bar(x + w/2, times_new, w, label=label_new, color=colors["p6_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)

    maxval = max(max(times_old), max(times_new))
    for bar in bars1:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + maxval * 0.015,
                f"{h:.2f}s", ha="center", va="bottom", fontsize=7)
    for bar in bars2:
        h = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2, h + maxval * 0.015,
                f"{h:.2f}s", ha="center", va="bottom", fontsize=7)

    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("执行时间 (s)")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, maxval * 1.3)
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


def plot_comparison_speedup(times_old, times_new, label_old, label_new, title, filename):
    sp_old = speedup(times_old)
    sp_new = speedup(times_new)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(threads, sp_old, "s-", color=colors["p5_speedup"], label=label_old, linewidth=1.8, markersize=6)
    ax.plot(threads, sp_new, "o-", color=colors["p6_speedup"], label=label_new, linewidth=2.2, markersize=7)

    for th, v in zip(threads, sp_old):
        ax.annotate(f"{v:.2f}x", (th, v), textcoords="offset points", xytext=(0, -12), ha="center", fontsize=7, color=colors["p5_speedup"])
    for th, v in zip(threads, sp_new):
        ax.annotate(f"{v:.2f}x", (th, v), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=7, color=colors["p6_speedup"])

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


def plot_all_phases_time(p1_t2, p5_t2, p6_t2, title, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    w = 0.25
    bars1 = ax.bar(x - w, p1_t2, w, label="阶段一 (Baseline)", color=colors["p1_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.7)
    bars2 = ax.bar(x,     p5_t2, w, label="阶段五 (Lock-Free)", color=colors["p5_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)
    bars3 = ax.bar(x + w, p6_t2, w, label="阶段六 (Greedy Early-Exit)", color=colors["p6_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)

    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("执行时间 (s)")
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper right")
    ax.grid(axis="y", alpha=0.25)
    ax.set_ylim(0, max(p1_t2) * 1.15)
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


def plot_all_phases_speedup(p1_t2, p5_t2, p6_t2, title, filename):
    base = p1_t2[0]
    sp1 = [base / v for v in p1_t2]
    sp5 = [base / v for v in p5_t2]
    sp6 = [base / v for v in p6_t2]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(threads, sp1, "^-", color=colors["p1_speedup"], label="阶段一", linewidth=1.5, markersize=6)
    ax.plot(threads, sp5, "s-", color=colors["p5_speedup"], label="阶段五", linewidth=1.8, markersize=6)
    ax.plot(threads, sp6, "o-", color=colors["p6_speedup"], label="阶段六 (Final)", linewidth=2.2, markersize=7)

    for th, v in zip(threads, sp6):
        ax.annotate(f"{v:.1f}x", (th, v), textcoords="offset points", xytext=(0, 10), ha="center", fontsize=8, fontweight="bold")

    ax.set_xlabel("线程数 (Thread Count)")
    ax.set_ylabel("加速比 (相对阶段一单线程)")
    ax.set_title(title)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log", base=2)
    ax.set_xticks(threads)
    ax.set_xticklabels(threads)
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.25, which="both")
    fig.savefig(os.path.join(OUTDIR, filename))
    plt.close(fig)


if __name__ == "__main__":
    plot_time(p6_times_t1, "阶段六执行时间 vs 线程数 (test.fasta, threshold=0.85)", "p6_exec_time_t1.png", "{:.3f}s")
    plot_speedup(p6_times_t1, "阶段六加速比曲线 (test.fasta)", "p6_speedup_t1.png")
    plot_efficiency(p6_times_t1, "阶段六并行效率 (test.fasta)", "p6_efficiency_t1.png")

    plot_time(p6_times_t2, "阶段六执行时间 vs 线程数 (test2.fasta, 四阈值平均)", "p6_exec_time_t2.png", "{:.3f}s")
    plot_speedup(p6_times_t2, "阶段六加速比曲线 (test2.fasta, 四阈值平均)", "p6_speedup_t2.png")
    plot_efficiency(p6_times_t2, "阶段六并行效率 (test2.fasta)", "p6_efficiency_t2.png")

    plot_comparison_time(p5_times_t2, p6_times_t2, "阶段五", "阶段六", "阶段五 vs 阶段六 执行时间 (test2.fasta)", "p6_vs_p5_exec_time.png")
    plot_comparison_speedup(p5_times_t2, p6_times_t2, "阶段五", "阶段六", "阶段五 vs 阶段六 加速比 (test2.fasta)", "p6_vs_p5_speedup.png")

    plot_all_phases_time(p1_times_t2, p5_times_t2, p6_times_t2, "全阶段执行时间对比 (test2.fasta)", "all_phases_exec_time_p6.png")
    plot_all_phases_speedup(p1_times_t2, p5_times_t2, p6_times_t2, "全阶段加速比对比 (相对阶段一单线程)", "all_phases_speedup_p6.png")

    print("Phase 6 charts generated successfully!")
