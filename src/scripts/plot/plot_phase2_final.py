"""
Generate Phase 2 performance charts for experiment report.
Data from remote run_all.ps1, 3 runs per thread count, all MD5 PASS.
Produces 6 charts:
  1. test.fasta  execution time
  2. test.fasta  speedup
  3. test.fasta  efficiency
  4. test2.fasta execution time
  5. test2.fasta speedup
  6. test2.fasta efficiency
  7. Phase 1 vs Phase 2 comparison (test2.fasta)
  8. Phase 1 vs Phase 2 speedup comparison (test2.fasta)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ── Phase 2 measured data ──────────────────────────────────────────
threads = [1, 2, 4, 8, 16, 32, 64]

# test.fasta (824 sequences)
p2_times_t1 = [0.073, 0.060, 0.040, 0.037, 0.030, 0.030, 0.040]

# test2.fasta (9697 sequences)
p2_times_t2 = [67.637, 34.963, 20.183, 12.537, 7.287, 4.683, 4.877]

# Phase 1 reference data (for comparison)
p1_times_t1 = [1.213, 0.653, 0.363, 0.213, 0.123, 0.107, 0.113]
p1_times_t2 = [184.910, 105.410, 56.573, 30.967, 18.547, 11.587, 9.713]

# ── Output directory ───────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
DEFAULT_OUTDIR = os.path.join(REPO_ROOT, "report_images")

outdir = sys.argv[1] if len(sys.argv) >= 2 else DEFAULT_OUTDIR
os.makedirs(outdir, exist_ok=True)

# ── Derived metrics ───────────────────────────────────────────────
def speedup(times):
    return [times[0] / t for t in times]

def efficiency(times):
    sp = speedup(times)
    return [sp[i] / threads[i] * 100 for i in range(len(threads))]

p2_sp_t1 = speedup(p2_times_t1)
p2_eff_t1 = efficiency(p2_times_t1)
p2_sp_t2 = speedup(p2_times_t2)
p2_eff_t2 = efficiency(p2_times_t2)

# ── Style ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans", "Arial"],
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

colors = {
    "time":       "#2d7dd2",
    "speedup":    "#27ae60",
    "ideal":      "#95a5a6",
    "efficiency": "#f39c12",
    "phase1":     "#bdc3c7",
    "phase2":     "#e74c3c",
}

x = np.arange(len(threads))

# ════════════════════════════════════════════════════════════════════
# Chart 1: test.fasta execution time
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, p2_times_t1, 0.55, color=colors["time"], edgecolor="white", linewidth=0.6)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + max(p2_times_t1) * 0.025,
            f"{h:.3f}s", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("执行时间 (s)")
ax.set_title("Phase 2 执行时间 vs 线程数 (test.fasta, threshold=0.85)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, max(p2_times_t1) * 1.25)
fig.savefig(os.path.join(outdir, "p2_exec_time_t1.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 2: test.fasta speedup
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
ax.plot(threads, p2_sp_t1, "o-", color=colors["speedup"], label="Phase 2 实测加速比",
        linewidth=2.2, markersize=7)
for th, sp in zip(threads, p2_sp_t1):
    ax.annotate(f"{sp:.2f}×", (th, sp), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8, color=colors["speedup"])
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("加速比 (Speedup)")
ax.set_title("Phase 2 加速比曲线 (test.fasta)")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc="upper left")
ax.grid(True, alpha=0.25, which="both")
fig.savefig(os.path.join(outdir, "p2_speedup_t1.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 3: test.fasta efficiency
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, p2_eff_t1, 0.55, color=colors["efficiency"], edgecolor="white", linewidth=0.6)
ax.axhline(y=100, color=colors["ideal"], linestyle="--", linewidth=1.2, label="理想效率 (100%)")
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 2,
            f"{h:.1f}%", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("并行效率 (%)")
ax.set_title("Phase 2 并行效率 (test.fasta)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, 120)
fig.savefig(os.path.join(outdir, "p2_efficiency_t1.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 4: test2.fasta execution time
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, p2_times_t2, 0.55, color=colors["time"], edgecolor="white", linewidth=0.6)
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + max(p2_times_t2) * 0.025,
            f"{h:.1f}s", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("执行时间 (s)")
ax.set_title("Phase 2 执行时间 vs 线程数 (test2.fasta, threshold=0.85)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, max(p2_times_t2) * 1.15)
fig.savefig(os.path.join(outdir, "p2_exec_time_t2.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 5: test2.fasta speedup
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.4)
ax.plot(threads, p2_sp_t2, "o-", color=colors["speedup"], label="Phase 2 实测加速比",
        linewidth=2.2, markersize=7)
for th, sp in zip(threads, p2_sp_t2):
    ax.annotate(f"{sp:.2f}×", (th, sp), textcoords="offset points", xytext=(0, 8),
                ha="center", fontsize=8, color=colors["speedup"])
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("加速比 (Speedup)")
ax.set_title("Phase 2 加速比曲线 (test2.fasta)")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc="upper left")
ax.grid(True, alpha=0.25, which="both")
fig.savefig(os.path.join(outdir, "p2_speedup_t2.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 6: test2.fasta efficiency
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(x, p2_eff_t2, 0.55, color=colors["efficiency"], edgecolor="white", linewidth=0.6)
ax.axhline(y=100, color=colors["ideal"], linestyle="--", linewidth=1.2, label="理想效率 (100%)")
for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 2,
            f"{h:.1f}%", ha="center", va="bottom", fontsize=8)
ax.set_xlabel("线程数 (Thread Count)")
ax.set_ylabel("并行效率 (%)")
ax.set_title("Phase 2 并行效率 (test2.fasta)")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.25)
ax.set_ylim(0, 120)
fig.savefig(os.path.join(outdir, "p2_efficiency_t2.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 7: Phase 1 vs Phase 2 execution time (test2.fasta)
# ════════════════════════════════════════════════════════════════════
fig, ax = plt.subplots(figsize=(10, 6))
w = 0.35
bars1 = ax.bar(x - w/2, p1_times_t2, w, color=colors["phase1"], edgecolor="white",
               linewidth=0.5, label="Phase 1")
bars2 = ax.bar(x + w/2, p2_times_t2, w, color=colors["phase2"], edgecolor="white",
               linewidth=0.5, label="Phase 2")

for bar in bars2:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, h + 1,
            f"{h:.1f}s", ha="center", va="bottom", fontsize=8, color="#c0392b", weight="bold")

ax.set_xlabel("线程数 (Thread Count)", fontsize=11, labelpad=10)
ax.set_ylabel("执行时间 (s)", fontsize=11, labelpad=10)
ax.set_title("Phase 1 vs Phase 2 执行时间对比 (test2.fasta)", fontsize=14, pad=15, weight="bold")
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc="upper right")
ax.grid(axis="y", alpha=0.3)

# Inset for high thread counts
axins = ax.inset_axes([0.45, 0.45, 0.4, 0.4])
axins.bar(x[-3:] - w/2, p1_times_t2[-3:], w, color=colors["phase1"])
axins.bar(x[-3:] + w/2, p2_times_t2[-3:], w, color=colors["phase2"])
axins.set_xticks(x[-3:])
axins.set_xticklabels(threads[-3:])
axins.set_title("16-64 线程放大", fontsize=10)
axins.grid(axis="y", alpha=0.3)
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.savefig(os.path.join(outdir, "p2_vs_p1_exec_time.png"))
plt.close(fig)

# ════════════════════════════════════════════════════════════════════
# Chart 8: Phase 1 vs Phase 2 speedup (test2.fasta)
# ════════════════════════════════════════════════════════════════════
p1_sp_t2 = speedup(p1_times_t2)

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(threads, threads, "--", color=colors["ideal"], label="理想加速比", linewidth=1.5)
ax.plot(threads, p1_sp_t2, "o--", color=colors["phase1"], label="Phase 1 加速比", linewidth=2)
ax.plot(threads, p2_sp_t2, "s-", color=colors["phase2"], label="Phase 2 加速比",
        linewidth=2.5, markersize=8)

for i, ps in enumerate(p2_sp_t2):
    ax.annotate(f"{ps:.1f}×", (threads[i], ps), textcoords="offset points", xytext=(0, -15),
                ha="center", fontsize=9, weight="bold", color=colors["phase2"])

ax.set_xlabel("线程数 (Thread Count)", fontsize=11, labelpad=10)
ax.set_ylabel("加速比 (Speedup)", fontsize=11, labelpad=10)
ax.set_title("Phase 1 vs Phase 2 加速比曲线 (test2.fasta)", fontsize=14, pad=15, weight="bold")
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc="upper left")
ax.grid(True, alpha=0.3, which="both")
fig.savefig(os.path.join(outdir, "p2_vs_p1_speedup.png"))
plt.close(fig)

print("All Phase 2 charts saved successfully!")
