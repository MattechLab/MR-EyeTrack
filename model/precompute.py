"""Precompute the k-space band used for training, one file per subject.

The trajectory is bit-identical across all 15 subjects (same .seq), so the
spoke directions are stored once and shared.

Band 160:320 keeps |k| <= 0.333 normalised -- structure down to ~6 mm, which
covers the globe (24 mm) and the lens (~9 mm). 818 MiB per subject.
"""
from pathlib import Path
import sys
import numpy as np
from dataset import load_band, load_labels, label_frame, CLASSES, k_axis

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
LO, HI = 160, 320

# Excluded: 2, 3, 5, 7 -- blink-event fractions of 0.42/0.26/0.22/0.21 are
# pupil-detection failure, not physiology (see model/DESIGN.md section 5).
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)


def main(subs=KEEP):
    OUT.mkdir(parents=True, exist_ok=True)
    for i, s in enumerate(subs):
        p = OUT / f'sub-{s:03d}_band{LO}-{HI}_ROI-PCA_8.npy'
        if p.exists():
            print(f'[{i+1}/{len(subs)}] sub-{s:03d} already done')
            continue
        y, d = load_band(s, LO, HI)
        np.save(p, y)
        if not (OUT / 'dirs.npy').exists():
            np.save(OUT / 'dirs.npy', d)
            np.save(OUT / 'kaxis.npy', k_axis(LO, HI))
        else:                                  # same sequence for every subject
            assert np.array_equal(d, np.load(OUT / 'dirs.npy'))
        L = label_frame(s)
        np.savez(OUT / f'sub-{s:03d}_labels.npz',
                 gaze=np.stack([load_labels(s, f'clean:{r}') for r in range(4)]),
                 **{c.replace('-', '_'): L[c] for c in CLASSES}, gap=L['gap'])
        print(f'[{i+1}/{len(subs)}] sub-{s:03d}  {y.shape}  {y.nbytes/2**20:.0f} MiB', flush=True)
    print('done ->', OUT)


if __name__ == '__main__':
    main([int(a) for a in sys.argv[1:]] or KEEP)
