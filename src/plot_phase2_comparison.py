import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

# Data
threads = [1, 2, 4, 8, 16, 32, 64]
phase1_times = [173.96, 89.85, 46.09, 24.47, 12.60, 8.52, 8.66]
phase2_times = [24.55, 12.72, 6.64, 3.55, 1.92, 1.23, 1.13]

outdir = sys.argv[1] if len(sys.argv) >= 2 else '.'

# Derived metrics
phase1_speedup = [phase1_times[0] / t for t in phase1_times]
phase2_speedup = [phase2_times[0] / t for t in phase2_times]
ideal_speedup = threads

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
    'phase1': '#bdc3c7',  # Grey for Phase 1
    'phase2': '#e74c3c',  # Vibrant red for Phase 2
    'ideal': '#34495e',   # Dark blue for ideal
}

# ============ Chart 1: Execution Time Comparison ============
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(threads))
w = 0.35

bars1 = ax.bar(x - w/2, phase1_times, w, color=colors['phase1'], edgecolor='white', linewidth=0.5, label='Phase 1 (精确计算)')
bars2 = ax.bar(x + w/2, phase2_times, w, color=colors['phase2'], edgecolor='white', linewidth=0.5, label='Phase 2 (TreeMinHash过滤)')

# Add text labels for Phase 2
for bar in bars2:
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2, height + 1,
            f'{height:.2f}s', ha='center', va='bottom', fontsize=8, color='#c0392b', weight='bold')

ax.set_xlabel('线程数 (Thread Count)', fontsize=11, labelpad=10)
ax.set_ylabel('执行时间 (s)', fontsize=11, labelpad=10)
ax.set_title('Phase 1 vs Phase 2 执行时间对比 (test2.fasta, N≈10000)', fontsize=14, pad=15, weight='bold')
ax.set_xticks(x)
ax.set_xticklabels(threads)
ax.legend(loc='upper right')
ax.grid(axis='y', alpha=0.3)
# ax.set_ylim(0, max(phase1_times) * 1.1)

# Inset plot for zoom on higher thread counts (16, 32, 64)
axins = ax.inset_axes([0.45, 0.45, 0.4, 0.4])
axins.bar(x[-3:] - w/2, phase1_times[-3:], w, color=colors['phase1'])
axins.bar(x[-3:] + w/2, phase2_times[-3:], w, color=colors['phase2'])
axins.set_xticks(x[-3:])
axins.set_xticklabels(threads[-3:])
axins.set_title('16-64 线程放大细节', fontsize=10)
axins.grid(axis='y', alpha=0.3)
ax.indicate_inset_zoom(axins, edgecolor="black")

fig.savefig(os.path.join(outdir, 'chart_exec_time_phase2.png'))
plt.close(fig)

# ============ Chart 2: Speedup Comparison ============
fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(threads, ideal_speedup, '--', color=colors['ideal'], label='理想加速比 (Ideal)', linewidth=1.5)
ax.plot(threads, phase1_speedup, 'o--', color=colors['phase1'], label='Phase 1 加速比', linewidth=2)
ax.plot(threads, phase2_speedup, 's-', color=colors['phase2'], label='Phase 2 加速比', linewidth=2.5, markersize=8)

for i, ps in enumerate(phase2_speedup):
    ax.annotate(f'{ps:.1f}x', (threads[i], ps), textcoords='offset points', xytext=(0, -15),
                ha='center', fontsize=9, weight='bold', color=colors['phase2'])

ax.set_xlabel('线程数 (Thread Count)', fontsize=11, labelpad=10)
ax.set_ylabel('加速比 (Speedup)', fontsize=11, labelpad=10)
ax.set_title('Phase 1 vs Phase 2 多线程加速比曲线 (test2.fasta)', fontsize=14, pad=15, weight='bold')
ax.set_xscale('log', base=2)
ax.set_yscale('log', base=2)
ax.set_xticks(threads)
ax.set_xticklabels(threads)
ax.legend(loc='upper left')
ax.grid(True, alpha=0.3, which="both")
fig.savefig(os.path.join(outdir, 'chart_speedup_phase2.png'))
plt.close(fig)

print("Phase 2 comparison charts saved successfully!")
