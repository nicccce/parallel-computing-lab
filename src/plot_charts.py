"""
Generate performance charts for the parallel computing experiment report.
Charts:
  1. Execution time vs. threads (baseline vs phase1)
  2. Speedup vs. threads (baseline vs phase1)
  3. Parallel efficiency vs. threads
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# ---- Data ----
threads = [1, 2, 4, 8, 16, 32, 64]

# Will be filled from command-line or defaults
baseline_times = None
phase1_times = None

def parse_args():
    """Parse benchmark output pasted as args, or use defaults."""
    global baseline_times, phase1_times
    if len(sys.argv) >= 3:
        baseline_times = list(map(float, sys.argv[1].split(',')))
        phase1_times   = list(map(float, sys.argv[2].split(',')))
    else:
        # Defaults from benchmark (will be overwritten)
        baseline_times = [1.15, 0.60, 0.35, 0.21, 0.15, 0.18, 0.20]
        phase1_times   = [1.36, 0.67, 0.38, 0.21, 0.14, 0.17, 0.19]

parse_args()

outdir = sys.argv[3] if len(sys.argv) >= 4 else '.'

# ---- Derived metrics ----
baseline_speedup = [baseline_times[0] / t for t in baseline_times]
phase1_speedup   = [phase1_times[0]   / t for t in phase1_times]
ideal_speedup    = threads

phase1_efficiency = [phase1_speedup[i] / threads[i] * 100 for i in range(len(threads))]
baseline_efficiency = [baseline_speedup[i] / threads[i] * 100 for i in range(len(threads))]

# ---- Style ----
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DejaVu Sans', 'Arial'],
    'axes.unicode_minus': False,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

colors = {
    'baseline': '#e74c3c',
    'phase1':   '#2ecc71',
    'ideal':    '#95a5a6',
}

# ============ Chart 1: Execution Time ============
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(threads))
w = 0.32

bars1 = ax.bar(x - w/2, baseline_times, w, label='Baseline (deepseek)', color=colors['baseline'], edgecolor='white', linewidth=0.5)
bars2 = ax.bar(x + w/2, phase1_times,   w, label='Phase 1 (optimized)', color=colors['phase1'],   edgecolor='white', linewidth=0.5)

for bar in bars1:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=7, color=colors['baseline'])
for bar in bars2:
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
            f'{bar.get_height():.2f}', ha='center', va='bottom', fontsize=7, color=colors['phase1'])

ax.set_xlabel('Thread Count')
ax.set_ylabel('Execution Time (s)')
ax.set_title('Execution Time vs. Thread Count (test.fasta, threshold=0.85)')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0, max(max(baseline_times), max(phase1_times)) * 1.25)
fig.savefig(os.path.join(outdir, 'chart_exec_time.png'))
plt.close(fig)

# ============ Chart 2: Speedup ============
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(threads, ideal_speedup, '--', color=colors['ideal'], label='Ideal Speedup', linewidth=1.5)
ax.plot(threads, baseline_speedup, 'o-', color=colors['baseline'], label='Baseline Speedup', linewidth=2, markersize=6)
ax.plot(threads, phase1_speedup,   's-', color=colors['phase1'],   label='Phase 1 Speedup',   linewidth=2, markersize=6)

for i, (bs, ps) in enumerate(zip(baseline_speedup, phase1_speedup)):
    ax.annotate(f'{bs:.1f}x', (threads[i], bs), textcoords='offset points', xytext=(0, 10),
                ha='center', fontsize=7, color=colors['baseline'])
    ax.annotate(f'{ps:.1f}x', (threads[i], ps), textcoords='offset points', xytext=(0, -15),
                ha='center', fontsize=7, color=colors['phase1'])

ax.set_xlabel('Thread Count')
ax.set_ylabel('Speedup')
ax.set_title('Speedup vs. Thread Count')
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig(os.path.join(outdir, 'chart_speedup.png'))
plt.close(fig)

# ============ Chart 3: Parallel Efficiency ============
fig, ax = plt.subplots(figsize=(8, 5))

ax.bar(np.arange(len(threads)) - w/2, baseline_efficiency, w,
       label='Baseline Efficiency', color=colors['baseline'], alpha=0.8, edgecolor='white')
ax.bar(np.arange(len(threads)) + w/2, phase1_efficiency, w,
       label='Phase 1 Efficiency',   color=colors['phase1'],   alpha=0.8, edgecolor='white')
ax.axhline(y=100, color=colors['ideal'], linestyle='--', linewidth=1, label='Ideal (100%)')

for i, (be, pe) in enumerate(zip(baseline_efficiency, phase1_efficiency)):
    ax.text(i - w/2, be + 1.5, f'{be:.0f}%', ha='center', va='bottom', fontsize=7, color=colors['baseline'])
    ax.text(i + w/2, pe + 1.5, f'{pe:.0f}%', ha='center', va='bottom', fontsize=7, color=colors['phase1'])

ax.set_xlabel('Thread Count')
ax.set_ylabel('Parallel Efficiency (%)')
ax.set_title('Parallel Efficiency vs. Thread Count')
ax.set_xticks(np.arange(len(threads)))
ax.set_xticklabels(threads)
ax.legend()
ax.grid(axis='y', alpha=0.3)
ax.set_ylim(0, 120)
fig.savefig(os.path.join(outdir, 'chart_efficiency.png'))
plt.close(fig)

print("Charts saved successfully!")
print(f"  - {os.path.join(outdir, 'chart_exec_time.png')}")
print(f"  - {os.path.join(outdir, 'chart_speedup.png')}")
print(f"  - {os.path.join(outdir, 'chart_efficiency.png')}")
