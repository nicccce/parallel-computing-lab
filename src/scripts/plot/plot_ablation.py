import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, "..", ".."))
OUTDIR = os.path.join(REPO_ROOT, "report_images")
os.makedirs(OUTDIR, exist_ok=True)

labels = ['Base (默认关闭)', '仅开启 AVX2 (Phase 4)', '仅开启倒排 (Phase 3)', '倒排 + AVX2 (全开)']
times = [20.912, 21.150, 20.895, 20.970]

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans", "Arial"],
    "axes.unicode_minus": False,
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "savefig.bbox": "tight",
})

fig, ax = plt.subplots(figsize=(8, 5))
x = np.arange(len(labels))
bars = ax.bar(x, times, width=0.5, color=['#2d7dd2', '#e74c3c', '#27ae60', '#f39c12'], edgecolor='white', linewidth=0.6)

for bar in bars:
    h = bar.get_height()
    ax.text(bar.get_x() + bar.get_width() / 2, h + 0.1, f"{h:.3f}s", ha="center", va="bottom", fontsize=10)

ax.set_ylabel("执行时间 (s)")
ax.set_title("微观消融实验对比 (4 Threads, test2.fasta 四阈值平均)")
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=15)
ax.set_ylim(20.0, 21.5)  # Zoom in to show the micro-differences
ax.grid(axis="y", alpha=0.3)

fig.savefig(os.path.join(OUTDIR, "ablation_study.png"))
plt.close(fig)
print("Ablation plot saved.")
