"""
Generate blink-comparison GIFs for each subject.

Per subject, two GIFs are saved under data/study/sub-XXX/recon/clean/:
  sub-XXX_axial_sl<n>_left_right.gif    — left (reg 2) vs right (reg 3)
  sub-XXX_sagittal_sl<n>_up_down.gif    — up   (reg 0) vs down  (reg 1)

Each GIF alternates between the two images at the specified slice so the
difference is visible as a blink comparison.
"""

import io
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

# ── Config ────────────────────────────────────────────────────────────────────
BASE_DIR  = "/home/debi/jaime/repos/MR-EyeTrack"
MASK_TYPE = "clean"
CLIP      = 0.7        # normalisation upper clip
DURATION  = 300        # ms per frame
DPI       = 120

# subject → (axial_slice, sagittal_slice)
SLICES = {
     1: (122, 147),
     2: (123, 146),
     3: (122,  94),
     4: (120,  97),
     5: (121,  94),
     6: (118,  95),
     7: (122, 144),
     8: (121,  88),
     9: (122, 147),
    10: (117, 146),
    11: (121, 146),
    12: (119, 144),
    13: (121, 152),
    14: (121,  95),
    15: (121, 146),
}

# ── Helpers ───────────────────────────────────────────────────────────────────
def norm(v):
    v = v.astype(float)
    return (v - v.min()) / (v.max() - v.min())

def norm_clip(v):
    return np.clip(norm(v), 0, CLIP)

def get_axial(vol, sl):
    return np.rot90(vol[:, :, sl], k=-1)

def to_sagittal(v):
    return np.transpose(np.rot90(v, k=-1, axes=(1, 2)), (1, 2, 0))

def fig_to_pil(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=DPI, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    return Image.open(buf).convert("RGBA")

def make_frame(slice_2d, title, label, cmap="gray", vmax=CLIP):
    fig, ax = plt.subplots(figsize=(4, 4), facecolor="black")
    ax.imshow(slice_2d, cmap=cmap, vmin=0, vmax=vmax, origin="lower")
    ax.set_title(label, color="white", fontsize=11, pad=4)
    ax.axis("off")
    fig.text(0.5, 0.01, title, ha="center", color="#aaaaaa", fontsize=8)
    img = fig_to_pil(fig)
    plt.close(fig)
    return img

def save_gif(frames, path):
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        loop=0,
        duration=DURATION,
        disposal=2,
    )
    print(f"  saved → {path}")

# ── Main loop ─────────────────────────────────────────────────────────────────
nii_tmpl = f"{BASE_DIR}/data/study/sub-{{s:03d}}/recon/{MASK_TYPE}/x/x_steva_regionidx_{{r}}_nIter_20_delta_1.000.nii.gz"
out_tmpl = f"{BASE_DIR}/data/study/sub-{{s:03d}}/recon/{MASK_TYPE}"

for subj, (ax_sl, sag_sl) in SLICES.items():
    print(f"\nsub-{subj:03d}  axial={ax_sl}  sagittal={sag_sl}")
    out_dir = out_tmpl.format(s=subj)

    # ── Axial: left (2) vs right (3) ─────────────────────────────────────────
    vol_l = norm_clip(np.abs(nib.load(nii_tmpl.format(s=subj, r=2)).get_fdata()))
    vol_r = norm_clip(np.abs(nib.load(nii_tmpl.format(s=subj, r=3)).get_fdata()))

    sl_l = get_axial(vol_l, ax_sl)
    sl_r = get_axial(vol_r, ax_sl)

    title_ax = f"sub-{subj:03d}  axial slice {ax_sl}  |  left (2) vs right (3)"
    frame_l = make_frame(sl_l, title_ax, "Left gaze (2)")
    frame_r = make_frame(sl_r, title_ax, "Right gaze (3)")

    gif_ax = f"{out_dir}/sub-{subj:03d}_axial_sl{ax_sl}_left_right.gif"
    save_gif([frame_l, frame_r], gif_ax)

    # ── Sagittal: up (0) vs down (1) ─────────────────────────────────────────
    vol_u = norm_clip(np.abs(nib.load(nii_tmpl.format(s=subj, r=0)).get_fdata()))
    vol_d = norm_clip(np.abs(nib.load(nii_tmpl.format(s=subj, r=1)).get_fdata()))

    sag_u = to_sagittal(vol_u)
    sag_d = to_sagittal(vol_d)

    sl_u = sag_u[:, :, sag_sl]
    sl_d = sag_d[:, :, sag_sl]

    title_sag = f"sub-{subj:03d}  sagittal slice {sag_sl}  |  up (0) vs down (1)"
    frame_u = make_frame(sl_u, title_sag, "Up gaze (0)")
    frame_d = make_frame(sl_d, title_sag, "Down gaze (1)")

    gif_sag = f"{out_dir}/sub-{subj:03d}_sagittal_sl{sag_sl}_up_down.gif"
    save_gif([frame_u, frame_d], gif_sag)

print("\nDone.")
