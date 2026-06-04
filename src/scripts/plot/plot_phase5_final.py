"""
Generate Phase 5 comprehensive tuning charts for experiment_report.md.
Data will be filled in from remote scripts/run/sweep_threads.ps1 results.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


threads = [1, 2, 4, 8, 16, 32, 64]

# ============================================================
# FILL IN from remote sweep_threads.ps1 results (all MD5 PASS)
# ============================================================

# test.fasta, threshold=0.85
p5_times_t1 = None  # Will be filled after benchmarking

# test2.fasta, average across thresholds 0.80/0.85/0.90/0.95
p5_times_t2 = None  # Will be filled after benchmarking

# Phase 4 data for comparison
p4_times_t1 = [0.080, 0.060, 0.040, 0.040, 0.030, 0.040, 0.040]
p4_times_t2 = [65.962, 34.812, 21.228, 12.845, 8.295, 5.620, 5.778]

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
    "time": "#2d7dd2",
    "speedup": "#27ae60",
    "ideal": "#95a5a6",
    "efficiency": "#f39c12",
    "p4_time": "#e74c3c",
    "p5_time": "#2d7dd2",
    "p4_speedup": "#e74c3c",
    "p5_speedup": "#27ae60",
    "p1_time": "#8e44ad",
    "p1_speedup": "#8e44ad",
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
    ax.plot(threads, sp, "o-", color=colors["speedup"], label="阶段五实测加速比", linewidth=2.2, markersize=7)
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
    """Side-by-side execution time comparison."""
    fig, ax = plt.subplots(figsize=(10, 5))
    w = 0.35
    bars1 = ax.bar(x - w/2, times_old, w, label=label_old, color=colors["p4_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)
    bars2 = ax.bar(x + w/2, times_new, w, label=label_new, color=colors["p5_time"],
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
    """Speedup comparison on log-log scale."""
    sp_old = speedup(times_old)
    sp_new = speedup(times_new)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(threads, sp_old, "s-", color=colors["p4_speedup"], label=label_old, linewidth=1.8, markersize=6)
    ax.plot(threads, sp_new, "o-", color=colors["p5_speedup"], label=label_new, linewidth=2.2, markersize=7)

    for th, v in zip(threads, sp_old):
        ax.annotate(f"{v:.2f}x", (th, v), textcoords="offset points", xytext=(0, -12), ha="center", fontsize=7, color=colors["p4_speedup"])
    for th, v in zip(threads, sp_new):
        ax.annotate(f"{v:.2f}x", (th, v), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=7, color=colors["p5_speedup"])

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


def plot_all_phases_time(p1_t2, p4_t2, p5_t2, title, filename):
    """All phases execution time comparison for test2.fasta."""
    fig, ax = plt.subplots(figsize=(10, 6))
    w = 0.25
    bars1 = ax.bar(x - w, p1_t2, w, label="阶段一 (Baseline)", color=colors["p1_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.7)
    bars2 = ax.bar(x,     p4_t2, w, label="阶段四", color=colors["p4_time"],
                   edgecolor="white", linewidth=0.6, alpha=0.85)
    bars3 = ax.bar(x + w, p5_t2, w, label="阶段五 (Final)", color=colors["p5_time"],
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


def plot_all_phases_speedup(p1_t2, p4_t2, p5_t2, title, filename):
    """All phases speedup comparison relative to Phase 1 single-thread."""
    base = p1_t2[0]
    sp1 = [base / v for v in p1_t2]
    sp4 = [base / v for v in p4_t2]
    sp5 = [base / v for v in p5_t2]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
    ax.plot(threads, sp1, "^-", color=colors["p1_speedup"], label="阶段一", linewidth=1.5, markersize=6)
    ax.plot(threads, sp4, "s-", color=colors["p4_speedup"], label="阶段四", linewidth=1.8, markersize=6)
    ax.plot(threads, sp5, "o-", color=colors["p5_speedup"], label="阶段五 (Final)", linewidth=2.2, markersize=7)

    for th, v in zip(threads, sp5):
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
    # Parse command-line data if provided:
    # Usage: python plot_phase5_final.py "t1_0,t1_1,...,t1_6" "t2_0,t2_1,...,t2_6"
    if len(sys.argv) >= 3:
        p5_times_t1 = [float(v) for v in sys.argv[1].split(",")]
        p5_times_t2 = [float(v) for v in sys.argv[2].split(",")]
    else:
        print("Usage: python plot_phase5_final.py <t1_csv> <t2_csv>")
        print("  e.g.: python plot_phase5_final.py '0.08,0.06,0.04,0.04,0.03,0.04,0.04' '64.0,33.5,20.0,12.0,8.0,5.5,5.6'")
        sys.exit(1)

    # Phase 5 standalone charts
    plot_time(p5_times_t1, "阶段五执行时间 vs 线程数 (test.fasta, threshold=0.85)",
              "p5_exec_time_t1.png", "{:.3f}s")
    plot_speedup(p5_times_t1, "阶段五加速比曲线 (test.fasta)", "p5_speedup_t1.png")
    plot_efficiency(p5_times_t1, "阶段五并行效率 (test.fasta)", "p5_efficiency_t1.png")

    plot_time(p5_times_t2, "阶段五执行时间 vs 线程数 (test2.fasta, 四阈值平均)",
              "p5_exec_time_t2.png", "{:.2f}s")
    plot_speedup(p5_times_t2, "阶段五加速比曲线 (test2.fasta, 四阈值平均)", "p5_speedup_t2.png")
    plot_efficiency(p5_times_t2, "阶段五并行效率 (test2.fasta)", "p5_efficiency_t2.png")

    # Phase 4 vs Phase 5 comparison
    plot_comparison_time(p4_times_t2, p5_times_t2,
                         "阶段四", "阶段五",
                         "阶段四 vs 阶段五 执行时间 (test2.fasta, 四阈值平均)",
                         "p5_vs_p4_exec_time.png")
    plot_comparison_speedup(p4_times_t2, p5_times_t2,
                            "阶段四", "阶段五",
                            "阶段四 vs 阶段五 加速比 (test2.fasta, 四阈值平均)",
                            "p5_vs_p4_speedup.png")

    # All phases comparison
    plot_all_phases_time(p1_times_t2, p4_times_t2, p5_times_t2,
                         "全阶段执行时间对比 (test2.fasta, 四阈值平均)",
                         "all_phases_exec_time.png")
    plot_all_phases_speedup(p1_times_t2, p4_times_t2, p5_times_t2,
                            "全阶段加速比对比 (test2.fasta, 相对阶段一单线程)",
                            "all_phases_speedup.png")

    print("All Phase 5 charts saved successfully!")
