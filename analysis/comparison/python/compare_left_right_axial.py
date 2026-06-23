"""
Compare left (region 2) vs right (region 3) — axial view.

Four panels: Left | Right | Difference | Toggle (blink)
  ↑ / ↓       move through slices
  ← / →       switch the toggle panel between Left and Right
"""

import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Slider

# ── Config ────────────────────────────────────────────────────────────────────
SUBJECT_NUM = 3
MASK_TYPE   = "clean"
BASE_DIR    = "/home/debi/jaime/repos/MR-EyeTrack"

CLIP    = 0.7   # upper clip after [0,1] normalisation
SL_INIT = 106   # initial axial slice

# ── Paths ─────────────────────────────────────────────────────────────────────
tmpl = f"{BASE_DIR}/data/study/sub-{SUBJECT_NUM:03d}/recon/{MASK_TYPE}/x/x_steva_regionidx_{{r}}_nIter_20_delta_1.000.nii.gz"
path_left  = tmpl.format(r=2)
path_right = tmpl.format(r=3)
print(f"left  path: {path_left}")
print(f"right path: {path_right}")

# ── Load ──────────────────────────────────────────────────────────────────────
vol_left  = np.abs(nib.load(path_left).get_fdata())
vol_right = np.abs(nib.load(path_right).get_fdata())

# ── Normalise ─────────────────────────────────────────────────────────────────
def norm(v):
    v = v.astype(float)
    return (v - v.min()) / (v.max() - v.min())

def norm_clip(v, clip=CLIP):
    return np.clip(norm(v), 0, clip)

img_left  = norm_clip(vol_left)
img_right = norm_clip(vol_right)

n_slices = img_left.shape[2]
SL_INIT  = max(0, min(SL_INIT, n_slices - 1))

diff_vol = norm(img_left) - norm(img_right)
diff_lim = np.abs(diff_vol).max() * 0.75

# ── Toggle state: 0 = show left, 1 = show right ───────────────────────────────
toggle_state = [0]
TOGGLE_LABELS = ["Toggle: Left (2)", "Toggle: Right (3)"]

# ── Figure ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(16, 5), facecolor="black")
fig.suptitle(
    f"sub-{SUBJECT_NUM:03d}  |  left (region 2) vs right (region 3)  |  axial\n"
    "↑↓ slices   ←→ toggle panel",
    color="white", fontsize=10
)

gs = gridspec.GridSpec(2, 4, figure=fig, height_ratios=[20, 1], hspace=0.08, wspace=0.05)

ax_left   = fig.add_subplot(gs[0, 0])
ax_right  = fig.add_subplot(gs[0, 1])
ax_diff   = fig.add_subplot(gs[0, 2])
ax_toggle = fig.add_subplot(gs[0, 3])

for ax, title in zip([ax_left, ax_right, ax_diff, ax_toggle],
                     ["Left (2)", "Right (3)", "Left − Right", TOGGLE_LABELS[0]]):
    ax.set_facecolor("black")
    ax.set_title(title, color="white", fontsize=9)
    ax.axis("off")

def get_axial(vol, sl):
    # rotate 90° counter-clockwise for display
    return np.rot90(vol[:, :, sl], k=-1)

def get_slice(sl):
    s_l    = get_axial(img_left,  sl)
    s_r    = get_axial(img_right, sl)
    s_diff = get_axial(diff_vol,  sl)
    s_tog  = s_l if toggle_state[0] == 0 else s_r
    return s_l, s_r, s_diff, s_tog

sl = SL_INIT
s_l, s_r, s_diff, s_tog = get_slice(sl)

im_left   = ax_left.imshow(s_l,    cmap="gray",   vmin=0, vmax=CLIP,     origin="lower")
im_right  = ax_right.imshow(s_r,   cmap="gray",   vmin=0, vmax=CLIP,     origin="lower")
im_diff   = ax_diff.imshow(s_diff, cmap="RdBu_r", vmin=-diff_lim, vmax=diff_lim, origin="lower")
im_toggle = ax_toggle.imshow(s_tog, cmap="gray",  vmin=0, vmax=CLIP,     origin="lower")

cbar_ax = fig.add_axes([0.51, 0.15, 0.01, 0.72])
fig.colorbar(im_diff, cax=cbar_ax)
cbar_ax.yaxis.set_tick_params(color="white", labelcolor="white")

# Slice slider
ax_slider = fig.add_subplot(gs[1, :])
slider = Slider(ax_slider, f"Axial slice  [0–{n_slices-1}]", 0, n_slices - 1,
                valinit=sl, valstep=1, color="#4a90d9")
ax_slider.set_facecolor("#1a1a1a")
slider.label.set_color("white")
slider.valtext.set_color("white")

def redraw():
    sl = int(slider.val)
    s_l, s_r, s_diff, s_tog = get_slice(sl)
    im_left.set_data(s_l)
    im_right.set_data(s_r)
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
