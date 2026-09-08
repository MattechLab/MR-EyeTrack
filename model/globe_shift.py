"""Did the eye actually rotate between gaze bins? A confound-free test.

The globe rotates about its own centre, so the lens (anterior) and the optic
nerve head (posterior) must move in **opposite** directions. Nothing else does
that: head translation, drift, gradient error and B0 shift all move the whole
orbit the same way. So the *differential* lens-minus-ONH displacement is immune
to every nuisance that has confounded the k-space analyses.

Predicted, from the stimulus geometry (visual_stimuli/mreyetrack_4points.py):

    left <-> right   13.75 deg  ->  lens ~2.4 mm, ONH ~2.9 mm opposite  (~5.3 mm apart)
    up   <-> down     7.78 deg  ->  lens ~1.4 mm, ONH ~1.6 mm opposite  (~3.0 mm apart)

Displacement is measured by sub-voxel cross-correlation of small boxes, so no
globe or lens segmentation is required.
"""
import sys
import numpy as np
import scipy.io as sio
from pathlib import Path

BASE = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
R_LENS, R_ONH = 10.0, 11.0          # mm from globe centre, anterior / posterior
BOX = 9                              # half-size of the correlation box, voxels


def steva(sub, r):
    p = BASE / f'sub-{sub:03d}/recon/clean/x/x_steva_regionidx_{r}_nIter_20_delta_1.000.mat'
    d = sio.loadmat(p)
    k = [x for x in d if not x.startswith('__')][0]
    return np.abs(np.asarray(d[k]).squeeze()).astype(np.float32)


def fit_globe(V, guess, search=6, radii=(10, 11, 12, 13)):
    """Globe centre = the sphere whose interior is most homogeneous and most
    distinct from the shell just outside it."""
    best = None
    g = np.array(guess)
    zz, yy, xx = np.mgrid[-16:17, -16:17, -16:17]
    d = np.sqrt(zz ** 2 + yy ** 2 + xx ** 2)
    for dz in range(-search, search + 1):
        for dy in range(-search, search + 1):
            for dx in range(-search, search + 1):
                c = g + [dz, dy, dx]
                sub = V[c[0]-16:c[0]+17, c[1]-16:c[1]+17, c[2]-16:c[2]+17]
                if sub.shape != d.shape:
                    continue
                for R in radii:
                    inside = sub[d <= R - 1]
                    shell = sub[(d > R + 1) & (d <= R + 4)]
                    if inside.size < 100 or shell.size < 100:
                        continue
                    score = abs(inside.mean() - shell.mean()) / (inside.std() + 1e-6)
                    if best is None or score > best[0]:
                        best = (score, tuple(c), R)
    return best


def shift(A, B, rad=5):
    """Sub-voxel displacement of B relative to A by normalised cross-correlation."""
    A = (A - A.mean()) / (A.std() + 1e-9)
    B = (B - B.mean()) / (B.std() + 1e-9)
    n = 2 * rad + 1
    cc = np.full((n, n, n), -np.inf)
    s = BOX
    for i, dz in enumerate(range(-rad, rad + 1)):
        for j, dy in enumerate(range(-rad, rad + 1)):
            for k, dx in enumerate(range(-rad, rad + 1)):
                b = B[s + dz - BOX + 1:s + dz + BOX, s + dy - BOX + 1:s + dy + BOX,
                      s + dx - BOX + 1:s + dx + BOX]
                a = A[1:2 * BOX, 1:2 * BOX, 1:2 * BOX]
                if b.shape != a.shape:
                    continue
                cc[i, j, k] = float((a * b).mean())
    p = np.unravel_index(np.argmax(cc), cc.shape)
    out = []
    for ax in range(3):
        i = p[ax]
        if 0 < i < n - 1:
            sl = [p[0], p[1], p[2]]
            sl[ax] = i - 1; ym = cc[tuple(sl)]
            sl[ax] = i + 1; yp = cc[tuple(sl)]
            y0 = cc[p]
            den = (ym - 2 * y0 + yp)
            off = 0.5 * (ym - yp) / den if abs(den) > 1e-12 else 0.0
        else:
            off = 0.0
        out.append(i - rad + off)
    return np.array(out)


def box(V, c):
    c = np.round(c).astype(int)
    return V[c[0]-BOX:c[0]+BOX+1, c[1]-BOX:c[1]+BOX+1, c[2]-BOX:c[2]+BOX+1]


def main():
    sub = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    V = [steva(sub, r) for r in range(4)]
    M = np.mean(V, 0)
    print(f'sub-{sub:03d}: axis0 = A-P (increasing = posterior), axis1 = L-R, axis2 = S-I\n')
    for eye, guess in (('left', (56, 95, 122)), ('right', (56, 152, 122))):
        sc, c, R = fit_globe(M, guess)
        c = np.array(c, float)
        print(f'{eye} globe: centre {c.astype(int)} radius {R} mm  (contrast score {sc:.2f})')
        lens = c + np.array([-R_LENS, 0, 0])          # anterior pole
        onh = c + np.array([+R_ONH, 0, 0])            # posterior pole
        ctrl = c + np.array([0, 0, 34])               # away from the orbit: should not move
        for nm, (a, b), axis, pred in (('left/right', (2, 3), 1, 2.4),
                                       ('up/down', (0, 1), 2, 1.4)):
            dl = shift(box(V[a], lens), box(V[b], lens))
            do = shift(box(V[a], onh), box(V[b], onh))
            dc = shift(box(V[a], ctrl), box(V[b], ctrl))
            axname = 'L-R' if axis == 1 else 'S-I'
            print(f'  {nm:11s} along {axname}:  lens {dl[axis]:+6.2f}  ONH {do[axis]:+6.2f}  '
                  f'control {dc[axis]:+6.2f}   differential {dl[axis]-do[axis]:+6.2f} mm '
                  f'(predicted {2*pred:+.1f})')
        print()
    print('The signature to look for: lens and ONH displacements with OPPOSITE sign,')
    print('a control box near zero, and the differential close to the prediction.')
    print('Head motion or drift would move lens, ONH and control the SAME way.')


if __name__ == '__main__':
    main()
