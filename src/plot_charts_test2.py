"""
Generate performance charts for the parallel Jaccard clustering on test2.fasta (Phase 1 optimized).
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

threads = [1, 2, 4, 8, 16, 32, 64]
phase1_times = [173.96, 89.85, 46.09, 24.47, 12.60, 8.52, 8.66]

outdir = sys.argv[1] if len(sys.argv) >= 2 else '.'

# Derived metrics
phase1_speedup = [phase1_times[0] / t for t in phase1_times]
ideal_speedup = threads
phase1_efficiency = [phase1_speedup[i] / threads[i] * 100 for i in range(len(threads))]

# Style
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DejaVu Sans', 'Arial'],
    'axes.unicode_minus': False,
    'figure.dpi': 150,
    'savefig.dpi': 150,
    'savefig.bbox': 'tight',
})

colors = {
    'phase1': '#3498db',  # Beautiful premium blue
    'ideal': '#bdc3c7',   # Elegant grey
    'accent': '#2ecc71',  # Emerald green for highlighting
}

# ============ Chart 1: Execution Time ============
fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(threads))
w = 0.5

bars = ax.bar(x, phase1_times, w, color=colors['phase1'], edgecolor='white', linewidth=0.5)

for bar in bars:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 3,
            f'{height:.2f}s', ha='center', va='bottom', fontsize=8, color='#2c3e50', weight='bold')

ax.set_xlabel('线程数 (Thread Count)', fontsize=10, labelpad=10)
ax.set_ylabel('执行时间 (s)', fontsize=10, labelpad=10)
ax.set_title('Phase 1 优化版执行时间与线程数关系 (test2.fasta, N≈10000)', fontsize=12, pad=15, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.grid(axis='y', alpha=0.2)
ax.set_ylim(0, max(phase1_times) * 1.15)
fig.savefig(os.path.join(outdir, 'chart_exec_time_test2.png'))
plt.close(fig)

# ============ Chart 2: Speedup ============
fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(threads, ideal_speedup, '--', color=colors['ideal'], label='理想加速比 (Ideal)', linewidth=1.5)
ax.plot(threads, phase1_speedup, 'o-', color=colors['phase1'], label='实际加速比 (Actual)', linewidth=2.5, markersize=8)

# Highlight some key points
for i, ps in enumerate(phase1_speedup):
    ax.annotate(f'{ps:.2f}x', (threads[i], ps), textcoords='offset points', xytext=(0, 10),
                ha='center', fontsize=8, weight='bold', color=colors['phase1'])

ax.set_xlabel('线程数 (Thread Count)', fontsize=10, labelpad=10)
ax.set_ylabel('加速比 (Speedup)', fontsize=10, labelpad=10)
ax.set_title('Phase 1 优化版加速比曲线 (test2.fasta)', fontsize=12, pad=15, weight='bold')
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.2, which="both")
fig.savefig(os.path.join(outdir, 'chart_speedup_test2.png'))
plt.close(fig)

# ============ Chart 3: Parallel Efficiency ============
fig, ax = plt.subplots(figsize=(8, 5))

bars_eff = ax.bar(x, phase1_efficiency, w, color='#2ecc71', alpha=0.85, edgecolor='white', linewidth=0.5)
ax.axhline(y=100, color=colors['ideal'], linestyle='--', linewidth=1.2, label='100% 理想效率')

for bar in bars_eff:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 2,
            f'{height:.1f}%', ha='center', va='bottom', fontsize=8, color='#27ae60', weight='bold')

ax.set_xlabel('线程数 (Thread Count)', fontsize=10, labelpad=10)
ax.set_ylabel('并行效率 (%)', fontsize=10, labelpad=10)
ax.set_title('Phase 1 优化版并行效率与线程数关系 (test2.fasta)', fontsize=12, pad=15, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc='lower left')
ax.grid(axis='y', alpha=0.2)
ax.set_ylim(0, 120)
fig.savefig(os.path.join(outdir, 'chart_efficiency_test2.png'))
plt.close(fig)

print("Charts for test2.fasta saved successfully!")
