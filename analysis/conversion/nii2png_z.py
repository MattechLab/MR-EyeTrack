import os
import sys
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons

# === Edit if needed ===
input_dir = "/home/debi/jaime/repos/MR-EyeTrack/results/bids/derivatives/reconstructions/sub-001/anat/"

# === Find stevaFull file (case-insensitive) ===
steva_fname = None
for fn in os.listdir(input_dir):
    if fn.endswith("stevaFull_T1w.nii.gz"):
        steva_fname = fn
        break

if steva_fname is None:
    print("ERROR: Could not find a file containing 'stevaFull' in the input directory.")
    print("Files found:")
    for fn in os.listdir(input_dir):
        print(" ", fn)
    sys.exit(1)

fpath = os.path.join(input_dir, steva_fname)
print("Using file:", fpath)

# === Load image (no reorientation, raw voxel space) ===
img = nib.load(fpath)
data = img.get_fdata()
nx, ny, nz = data.shape
print(f"Loaded {steva_fname} with shape {data.shape}")

# === Initial parameters ===
z0 = nz // 2
rotation_k = 0  # np.rot90(..., k) where k = 0,1,2,3 for 0°,90°,180°,270°

# === Figure and axes ===
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(left=0.15, bottom=0.25, right=0.9, top=0.92)

# Display the initial axial slice (we'll rotate later using k)
def get_rotated_axial(z_idx, k):
    slc = data[:, :, z_idx]
    if k != 0:
        slc = np.rot90(slc, k=k)
    return slc

im = ax.imshow(get_rotated_axial(z0, rotation_k), cmap="gray", origin="lower", interpolation="nearest")
ax.set_title(f"Axial (z = {z0})")
ax.axis("off")

# === Slider for Z ===
ax_slider = plt.axes([0.15, 0.12, 0.65, 0.03])  # [left, bottom, width, height]
z_slider = Slider(ax_slider, "Z", 0, nz - 1, valinit=z0, valstep=1, valfmt='%0.0f')

def on_slider(val):
    z_idx = int(z_slider.val)
    im.set_data(get_rotated_axial(z_idx, rotation_k))
    ax.set_title(f"Axial (z = {z_idx})")
    fig.canvas.draw_idle()

z_slider.on_changed(on_slider)

# === Radio buttons for rotation (0, 90, 180, 270) ===
ax_radio = plt.axes([0.82, 0.35, 0.12, 0.2], facecolor='lightgoldenrodyellow')
radio = RadioButtons(ax_radio, ("0°", "90°", "180°", "270°"), active=0)

def on_radio(label):
    global rotation_k
    if label == "0°":
        rotation_k = 0
    elif label == "90°":
        rotation_k = 1
    elif label == "180°":
        rotation_k = 2
    elif label == "270°":
        rotation_k = 3
    # update image using current slider z
    z_idx = int(z_slider.val)
    im.set_data(get_rotated_axial(z_idx, rotation_k))
    fig.canvas.draw_idle()

radio.on_clicked(on_radio)

# === Instructions text ===
instr = (
    "Use the slider to change Z (axial). \n"
    "Pick rotation to match Mango. \n"
    "When you find the matching slice, tell me the Z value shown in the title."
)
fig.text(0.15, 0.02, instr, fontsize=9)

print("GUI ready: use the slider to change z, and the radio buttons to set rotation (0/90/180/270°).")
plt.show()
