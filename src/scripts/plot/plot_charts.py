"""
Generate performance charts for the current optimized implementation on
test.fasta. Data points are the averages of three remote run_all.ps1 runs on
2026-06-04, all with MD5 PASS.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


threads = [1, 2, 4, 8, 16, 32, 64]
times = [1.213, 0.653, 0.363, 0.213, 0.123, 0.107, 0.113]

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "report_images")

outdir = sys.argv[1] if len(sys.argv) >= 2 else DEFAULT_OUTDIR
os.makedirs(outdir, exist_ok=True)

speedup = [times[0] / t for t in times]
efficiency = [speedup[i] / threads[i] * 100 for i in range(len(threads))]

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

# Execution time
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, times, 0.55, color=colors["time"], edgecolor="white", linewidth=0.6)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + max(times) * 0.025,
            f"{h:.3f}s", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("Thread Count")
ax.set_ylabel("Execution Time (s)")
ax.set_title("Execution Time vs. Thread Count (test.fasta, threshold=0.85)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, max(times) * 1.18)
fig.savefig(os.path.join(outdir, "chart_exec_time.png"))
plt.close(fig)

# Speedup
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(threads, threads, "--", color=colors["ideal"], label="Ideal Speedup", linewidth=1.4)
ax.plot(threads, speedup, "o-", color=colors["speedup"], label="Measured Speedup",
        linewidth=2.2, markersize=7)
for th, sp in zip(threads, speedup):
    ax.annotate(f"{sp:.2f}x", (th, sp), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8, color=colors["speedup"])
ax.set_xlabel("Thread Count")
ax.set_ylabel("Speedup")
ax.set_title("Speedup vs. Thread Count (test.fasta)")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc="upper left")
ax.grid(True, alpha=0.25, which="both")
fig.savefig(os.path.join(outdir, "chart_speedup.png"))
plt.close(fig)

# Efficiency
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, efficiency, 0.55, color=colors["efficiency"], edgecolor="white",
              linewidth=0.6)
ax.axhline(y=100, color=colors["ideal"], linestyle="--", linewidth=1.2,
           label="Ideal (100%)")
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 2,
            f"{h:.1f}%", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("Thread Count")
ax.set_ylabel("Parallel Efficiency (%)")
ax.set_title("Parallel Efficiency vs. Thread Count (test.fasta)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, 120)
fig.savefig(os.path.join(outdir, "chart_efficiency.png"))
plt.close(fig)

print("Charts for test.fasta saved successfully.")
