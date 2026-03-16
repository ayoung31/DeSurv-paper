"""
DeSurv hex sticker.
Portrait hexagon (pointy-top): 2 in × 2.309 in standard R hex sticker dimensions.
Rendered at 600 DPI for print quality.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from matplotlib.collections import PatchCollection
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patheffects as pe

# ── Palette ───────────────────────────────────────────────────────────────────
NAVY   = "#1B355E"
TEAL   = "#007A7C"
ORANGE = "#E07B1C"
WHITE  = "#FFFFFF"
LTBLUE = "#5BC8C8"
DIMTXT = "#AABBCC"

# ── Canvas ────────────────────────────────────────────────────────────────────
W, H = 4.0, 4.618
DPI  = 600

fig = plt.figure(figsize=(W, H), dpi=DPI)
ax  = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.set_aspect("equal")
ax.axis("off")
fig.patch.set_alpha(0)


# ── Hex geometry (pointy-top) ─────────────────────────────────────────────────
cx, cy = W / 2, H / 2
r = (W / np.sqrt(3)) * 0.960
angles = np.radians([30, 90, 150, 210, 270, 330])
hx = cx + r * np.cos(angles)
hy = cy + r * np.sin(angles)

_verts = list(zip(hx, hy)) + [(hx[0], hy[0])]
_codes = [Path.MOVETO] + [Path.LINETO] * 5 + [Path.CLOSEPOLY]
HEX = Path(_verts, _codes)


def hex_clip():
    return PathPatch(HEX, transform=ax.transData)


# ── 1. Navy background ────────────────────────────────────────────────────────
ax.add_patch(PathPatch(HEX, facecolor=NAVY, edgecolor="none", zorder=0))


# ── 2. Heatmap grid (more pronounced: higher alpha, stronger clusters) ────────
np.random.seed(42)
N_COLS, N_ROWS = 30, 11
CELL_W = W / N_COLS
CELL_H = (H * 0.32) / N_ROWS
HM_BOTTOM = H * 0.615

expr = np.zeros((N_ROWS, N_COLS))
# Strong cluster blocks + low background noise
expr[:4,   :10] += np.random.gamma(4.5, 1.0, (4, 10))   # cluster 1
expr[4:8,  10:20] += np.random.gamma(4.5, 1.0, (4, 10)) # cluster 2
expr[8:,   20:]   += np.random.gamma(4.5, 1.0, (N_ROWS - 8, N_COLS - 20))  # cluster 3
expr += np.random.gamma(0.3, 0.2, expr.shape)            # background noise
expr = (expr - expr.min()) / (expr.max() - expr.min())

cmap = LinearSegmentedColormap.from_list(
    "desurv_hm", ["#0A1628", "#0B4050", "#007A7C", "#D4691A"], N=256
)

rects, vals = [], []
for i in range(N_ROWS):
    for j in range(N_COLS):
        x = HM_BOTTOM + j * CELL_W   # intentional: using HM_BOTTOM as left offset = 0
        x = j * CELL_W
        y = HM_BOTTOM + i * CELL_H
        rects.append(mpatches.Rectangle((x, y), CELL_W * 0.92, CELL_H * 0.88))
        vals.append(expr[i, j])

pc = PatchCollection(rects, cmap=cmap, alpha=0.82, zorder=1, linewidth=0)
pc.set_array(np.array(vals))
pc.set_clim(0, 1)
pc.set_clip_path(hex_clip())
ax.add_collection(pc)


# ── 3. KM survival curves ─────────────────────────────────────────────────────
T_MAX   = 60
N_STEPS = 40
t = np.linspace(0, T_MAX, N_STEPS)

s_low  = np.exp(-0.016 * t)
s_high = np.exp(-0.052 * t)

KM_X0, KM_X1 = W * 0.10, W * 0.90
KM_Y0, KM_Y1 = H * 0.175, H * 0.455   # raised floor so URL sits below curves


def t2x(tv): return KM_X0 + (tv / T_MAX) * (KM_X1 - KM_X0)
def s2y(sv): return KM_Y0 + sv * (KM_Y1 - KM_Y0)


def make_steps(tv, sv):
    xs, ys = [], []
    for i in range(len(tv) - 1):
        xs += [t2x(tv[i]), t2x(tv[i + 1])]
        ys += [s2y(sv[i]), s2y(sv[i])]
    xs.append(t2x(tv[-1]))
    ys.append(s2y(sv[-1]))
    return xs, ys


xs_l, ys_l = make_steps(t, s_low)
xs_h, ys_h = make_steps(t, s_high)

ax.fill_between(
    [t2x(ti) for ti in t],
    [s2y(si) for si in s_high],
    [s2y(si) for si in s_low],
    color=WHITE, alpha=0.07, zorder=2, clip_path=hex_clip()
)

lw = 3.8
line_l, = ax.plot(xs_l, ys_l, color=LTBLUE, lw=lw,
                  solid_capstyle="round", solid_joinstyle="round", zorder=4)
line_h, = ax.plot(xs_h, ys_h, color=ORANGE, lw=lw,
                  solid_capstyle="round", solid_joinstyle="round", zorder=4)
line_l.set_clip_path(hex_clip())
line_h.set_clip_path(hex_clip())

# ── 4. Package name — measure then true-center ────────────────────────────────
# First pass: draw invisible proxy texts to measure pixel widths
FONT_SIZE = 46
name_y = H * 0.590

_m_de   = ax.text(0, 0, "De",   fontsize=FONT_SIZE, fontweight="bold",
                  fontfamily="DejaVu Sans", alpha=0, zorder=-1)
_m_surv = ax.text(0, 0, "Surv", fontsize=FONT_SIZE, fontweight="bold",
                  fontfamily="DejaVu Sans", alpha=0, zorder=-1)

fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv = ax.transData.inverted()

def px_width(txt_obj):
    bb = txt_obj.get_window_extent(renderer)
    return (inv.transform((bb.width, 0)) - inv.transform((0, 0)))[0]

GAP = W * 0.012   # small gap between De and Surv
w_de   = px_width(_m_de)
w_surv = px_width(_m_surv)
_m_de.remove()
_m_surv.remove()

total_w  = w_de + GAP + w_surv
de_cx    = cx - total_w / 2 + w_de / 2        # center of "De"
surv_cx  = cx + total_w / 2 - w_surv / 2      # center of "Surv"

t_de = ax.text(de_cx, name_y, "De",
               fontsize=FONT_SIZE, fontweight="bold", color=WHITE,
               ha="center", va="center", zorder=5,
               fontfamily="DejaVu Sans",
               path_effects=[pe.withStroke(linewidth=3, foreground=NAVY)])
t_de.set_clip_path(hex_clip())
t_de.set_clip_on(True)

t_surv = ax.text(surv_cx, name_y, "Surv",
                 fontsize=FONT_SIZE, fontweight="bold", color=ORANGE,
                 ha="center", va="center", zorder=5,
                 fontfamily="DejaVu Sans",
                 path_effects=[pe.withStroke(linewidth=3, foreground=NAVY)])
t_surv.set_clip_path(hex_clip())
t_surv.set_clip_on(True)


# ── 5. Tagline ────────────────────────────────────────────────────────────────
ax.text(W / 2, H * 0.500, "survival-guided NMF",
        fontsize=15.5, fontstyle="italic", color=LTBLUE,
        ha="center", va="center", zorder=5,
        fontfamily="DejaVu Sans",
        path_effects=[pe.withStroke(linewidth=2, foreground=NAVY)])


# ── 6. URL — backed by a translucent rect so it reads over the curves ─────────
URL_Y  = H * 0.143
URL_BH = H * 0.055
URL_BW = W * 0.62
url_bg = mpatches.FancyBboxPatch(
    (cx - URL_BW / 2, URL_Y - URL_BH / 2), URL_BW, URL_BH,
    boxstyle="round,pad=0.01",
    facecolor=NAVY, edgecolor="none", alpha=0.72, zorder=5
)
url_bg.set_clip_path(hex_clip())
ax.add_patch(url_bg)

t_url = ax.text(W / 2, URL_Y, "ayoung31/DeSurv",
                fontsize=11, color=DIMTXT,
                ha="center", va="center", zorder=6,
                fontfamily="DejaVu Sans")
t_url.set_clip_path(hex_clip())
t_url.set_clip_on(True)


# ── 7. Hex border ─────────────────────────────────────────────────────────────
ax.add_patch(PathPatch(HEX, facecolor="none",
                       edgecolor=TEAL, linewidth=18, zorder=7))
ax.add_patch(PathPatch(HEX, facecolor="none",
                       edgecolor=WHITE, linewidth=3, alpha=0.25, zorder=8))


# ── Save ──────────────────────────────────────────────────────────────────────
out = r"C:\Users\amyou\Documents\GitHub\DeSurv-paper\DeSurv_hex.png"
fig.savefig(out, dpi=DPI, bbox_inches="tight", pad_inches=0, transparent=True)
print(f"Saved: {out}")
plt.close(fig)
