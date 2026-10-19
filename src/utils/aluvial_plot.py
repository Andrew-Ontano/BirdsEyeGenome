import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.cm as cm

# ==========================================
# 1. Simulate Dummy Data (Replace with yours)
# ==========================================
np.random.seed(42)
num_windows = 32000
levels = list(range(2, 28))  # 26 levels (k=2 to k=27)

# Matrix shape: (26 levels, 32,000 windows)
# Rows = Levels (k=2 to 27), Columns = Genomic Position
data = np.zeros((len(levels), num_windows), dtype=int)
for idx, k in enumerate(levels):
    # Generating dummy cluster assignments (0 to k-1)
    data[idx, :] = np.random.randint(0, k, size=num_windows)

# ==========================================
# 2. Build a Discrete Qualitative Colormap
# ==========================================
max_k = 27
# Combine qualitative colormaps to get 27 distinct colors
colors_tab20 = [cm.get_cmap('tab20')(i) for i in range(20)]
colors_set3  = [cm.get_cmap('Set3')(i) for i in range(7)]
custom_cmap = ListedColormap(colors_tab20 + colors_set3)

# ==========================================
# 3. Plot Stacked Track Matrix
# ==========================================
fig, ax = plt.subplots(figsize=(16, 8))

im = ax.imshow(
    data,
    aspect='auto',            # Stretches y-axis to fit x-axis aspect ratio
    cmap=custom_cmap,
    interpolation='nearest',   # Prevents pixel blending between neighboring windows
    origin='lower',            # Lowest level (k=2) at the bottom
    vmin=0,
    vmax=max_k - 1,
    rasterized=True            # Critical: rasterizes heatmap pixels when saving vector PDFs
)

# ==========================================
# 4. Axes & Labels
# ==========================================
ax.set_title("Hierarchical Domain Assignments (32,000 Windows across k=2 to 27)", fontsize=14, pad=12)
ax.set_xlabel("Genomic Window Index (1 to 32,000)", fontsize=12)
ax.set_ylabel("Cluster Resolution Level", fontsize=12)

# Align Y-axis ticks with level labels
ax.set_yticks(range(len(levels)))
ax.set_yticklabels([f"k={k}" for k in levels], fontsize=8)

# Add colorbar for cluster IDs
cbar = fig.colorbar(im, ax=ax, pad=0.015, aspect=30)
cbar.set_label("Cluster ID Assignment", fontsize=11)

plt.tight_layout()

# Save as high-res PNG or PDF (rasterized=True prevents huge PDF file sizes)
plt.savefig("genomic_hierarchy_tracks.png", dpi=300, bbox_inches='tight')
plt.show()