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
    if fn.lower().endswith("stevafull_t1w.nii.gz"):
        steva_fname = fn
        break

if steva_fname is None:
    print("ERROR: Could not find a file containing 'stevaFull' in the input directory.")
    sys.exit(1)

fpath = os.path.join(input_dir, steva_fname)
print("Using file:", fpath)

# === Load image (raw voxel space) ===
img = nib.load(fpath)
data = img.get_fdata()
nx, ny, nz = data.shape
print(f"Loaded {steva_fname} with shape {data.shape}")

# === Initial parameters for X-axis (sagittal) ===
x0 = nx // 2
rotation_k = 1  # 90° default for sagittal Mango orientation

# === Figure and axes ===
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(left=0.15, bottom=0.25, right=0.9, top=0.92)

# Function to extract sagittal slice and apply rotation/flip
def get_sagittal_slice(x_idx, k):
    slc = data[x_idx, :, :]
    slc = np.rot90(slc, k)   # rotate
    slc = np.fliplr(slc)      # horizontal flip to match Mango left-right
    return slc

# Display initial sagittal slice
im = ax.imshow(get_sagittal_slice(x0, rotation_k), cmap="gray", origin="lower", interpolation="nearest")
ax.set_title(f"Sagittal (x = {x0})")
ax.axis("off")

# === Slider for X ===
ax_slider = plt.axes([0.15, 0.12, 0.65, 0.03])
x_slider = Slider(ax_slider, "X", 0, nx - 1, valinit=x0, valstep=1, valfmt='%0.0f')

def on_slider(val):
    x_idx = int(x_slider.val)
    im.set_data(get_sagittal_slice(x_idx, rotation_k))
    ax.set_title(f"Sagittal (x = {x_idx})")
    fig.canvas.draw_idle()

x_slider.on_changed(on_slider)

# === Radio buttons for rotation ===
ax_radio = plt.axes([0.82, 0.35, 0.12, 0.2], facecolor='lightgoldenrodyellow')
radio = RadioButtons(ax_radio, ("0°", "90°", "180°", "270°"), active=1)

def on_radio(label):
    global rotation_k
    rotation_k = {"0°":0, "90°":1, "180°":2, "270°":3}[label]
    x_idx = int(x_slider.val)
    im.set_data(get_sagittal_slice(x_idx, rotation_k))
    fig.canvas.draw_idle()

radio.on_clicked(on_radio)

# === Instructions text ===
instr = (
    "Use the slider to change X (sagittal). \n"
    "Pick rotation to match Mango orientation. \n"
    "When you find the matching slice, tell me the X value shown in the title."
)
fig.text(0.15, 0.02, instr, fontsize=9)

print("GUI ready: use the slider to change x, and the radio buttons to set rotation (0/90/180/270°).")
plt.show()
