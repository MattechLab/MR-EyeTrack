"""
Compare up (region 0) vs down (region 1) — sagittal view.

Slices along x (dim 0), same pattern as axial uses z (dim 2).
rot90 k=-1 on each (AP, SI) slice → SI vertical, anterior left.

Four panels: Up | Down | Difference | Toggle (blink)
  ↑ / ↓   move through sagittal slices (left-to-right across brain)
  ← / →   switch the toggle panel between Up and Down
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider

# ── Config ────────────────────────────────────────────────────────────────────
SUBJECT_NUM = 15
MASK_TYPE   = "clean"
BASE_DIR    = "/home/debi/jaime/repos/MR-EyeTrack"

CLIP    = 0.7   # upper clip after [0,1] normalisation
SL_INIT = 120   # initial sagittal slice (mid-brain left-right)

# ── Paths ─────────────────────────────────────────────────────────────────────
tmpl = f"{BASE_DIR}/data/study/sub-{SUBJECT_NUM:03d}/recon/{MASK_TYPE}/x/x_steva_regionidx_{{r}}_nIter_20_delta_1.000.nii.gz"
path_up   = tmpl.format(r=0)
path_down = tmpl.format(r=1)
print(f"up   path: {path_up}")
print(f"down path: {path_down}")

# ── Load ──────────────────────────────────────────────────────────────────────
vol_up   = np.abs(nib.load(path_up).get_fdata())
vol_down = np.abs(nib.load(path_down).get_fdata())

# ── Normalise ─────────────────────────────────────────────────────────────────
def norm(v):
    v = v.astype(float)
    return (v - v.min()) / (v.max() - v.min())

def norm_clip(v, clip=CLIP):
    return np.clip(norm(v), 0, clip)

img_up   = norm_clip(vol_up)
img_down = norm_clip(vol_down)

# ── Sagittal permutation ───────────────────────────────────────────────────────
# Data is RAS: dim0=R, dim1=A, dim2=S.
# Sagittal: slice along dim0 (x/R).  Each vol[i,:,:] is shape (Ny=AP, Nz=SI).
# rot90 k=-1 per-slice: SI vertical (inferior at bottom), anterior at left.
# np.rot90(v, k=-1, axes=(1,2)) applies this to the whole volume → (Nx, Nz, Ny).
# Transpose to (Nz, Ny, Nx) so sagittal slices are along the last axis.
def to_sagittal(v):
    return np.transpose(np.rot90(v, k=-1, axes=(1, 2)), (1, 2, 0))

sag_up   = to_sagittal(img_up)
sag_down = to_sagittal(img_down)

n_slices = sag_up.shape[2]   # = Nx = 240
SL_INIT  = max(0, min(SL_INIT, n_slices - 1))

diff_vol = norm(sag_up) - norm(sag_down)
diff_lim = np.abs(diff_vol).max() * 0.75

# ── Toggle state: 0 = show up, 1 = show down ─────────────────────────────────
toggle_state  = [0]
TOGGLE_LABELS = ["Toggle: Up (0)", "Toggle: Down (1)"]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 5), facecolor="black")
fig.suptitle(
    f"sub-{SUBJECT_NUM:03d}  |  up (region 0) vs down (region 1)  |  sagittal\n"
    "↑↓ slices   ←→ toggle panel",
    color="white", fontsize=10
)

gs = gridspec.GridSpec(2, 4, figure=fig, height_ratios=[20, 1], hspace=0.08, wspace=0.05)

ax_up     = fig.add_subplot(gs[0, 0])
ax_down   = fig.add_subplot(gs[0, 1])
ax_diff   = fig.add_subplot(gs[0, 2])
ax_toggle = fig.add_subplot(gs[0, 3])

for ax, title in zip([ax_up, ax_down, ax_diff, ax_toggle],
                     ["Up (0)", "Down (1)", "Up − Down", TOGGLE_LABELS[0]]):
    ax.set_facecolor("black")
    ax.set_title(title, color="white", fontsize=9)
    ax.axis("off")

def get_slice(sl):
    s_u    = sag_up[:, :, sl]
    s_d    = sag_down[:, :, sl]
    s_diff = diff_vol[:, :, sl]
    s_tog  = s_u if toggle_state[0] == 0 else s_d
    return s_u, s_d, s_diff, s_tog

sl = SL_INIT
s_u, s_d, s_diff, s_tog = get_slice(sl)

im_up     = ax_up.imshow(s_u,     cmap="gray",   vmin=0, vmax=CLIP,              origin="lower")
im_down   = ax_down.imshow(s_d,   cmap="gray",   vmin=0, vmax=CLIP,              origin="lower")
im_diff   = ax_diff.imshow(s_diff, cmap="RdBu_r", vmin=-diff_lim, vmax=diff_lim, origin="lower")
im_toggle = ax_toggle.imshow(s_tog, cmap="gray",  vmin=0, vmax=CLIP,             origin="lower")

cbar_ax = fig.add_axes([0.51, 0.15, 0.01, 0.72])
fig.colorbar(im_diff, cax=cbar_ax)
cbar_ax.yaxis.set_tick_params(color="white", labelcolor="white")

# Slice slider
ax_slider = fig.add_subplot(gs[1, :])
slider = Slider(ax_slider, f"Sagittal slice  [0–{n_slices-1}]", 0, n_slices - 1,
                valinit=sl, valstep=1, color="#4a90d9")
ax_slider.set_facecolor("#1a1a1a")
slider.label.set_color("white")
slider.valtext.set_color("white")

def redraw():
    sl = int(slider.val)
    s_u, s_d, s_diff, s_tog = get_slice(sl)
    im_up.set_data(s_u)
    im_down.set_data(s_d)
    im_diff.set_data(s_diff)
    im_toggle.set_data(s_tog)
    ax_toggle.set_title(TOGGLE_LABELS[toggle_state[0]], color="white", fontsize=9)
    fig.canvas.draw_idle()

slider.on_changed(lambda val: redraw())

def on_key(event):
    sl = int(slider.val)
    if event.key == "up":
        slider.set_val(min(sl + 1, n_slices - 1))
    elif event.key == "down":
        slider.set_val(max(sl - 1, 0))
    elif event.key in ("left", "right"):
        toggle_state[0] = 1 - toggle_state[0]
        redraw()

fig.canvas.mpl_connect("key_press_event", on_key)

plt.show()
