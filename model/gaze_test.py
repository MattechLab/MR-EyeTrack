"""Definitive test: is there per-readout gaze separation in the k-space?

d'^2 against a null built from the same statistic with the labels circularly
shifted. The shift preserves each class's block structure and duty cycle, so
the null absorbs everything except the actual time-locking of labels to data --
which is the only thing a real effect can be.
"""
import sys
import numpy as np
from pathlib import Path
from dprime import prep, dprime, OUT, NR

SHIFTS = (7000, 13000, 20000, 27000, 34000, 41000, 48000, 55000)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
CONTRASTS = {'up/down': (0, 1), 'left/right': (2, 3), 'horiz/vert': None}


def test(sub, kpc=32):
    Z = prep(sub, kpc=kpc)
    g = np.load(OUT / f'sub-{sub:03d}_labels.npz')['gaze']
    rows = {}
    for name, pair in CONTRASTS.items():
        if pair is None:
            a, b = g[0] | g[1], g[2] | g[3]
        else:
            a, b = g[pair[0]], g[pair[1]]
        obs, nuse = dprime(Z, a, b)
        nul = np.array([dprime(Z, a, b, shift=s)[0] for s in SHIFTS])
        z = (obs - nul.mean()) / (nul.std(ddof=1) + 1e-12)
        rows[name] = (obs, nul.mean(), nul.std(ddof=1), z, nuse)
    return rows


if __name__ == '__main__':
    kpc = 32
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    print(f"per-readout d'^2, gaze contrasts, {kpc} PCs, 8-shift null\n")
    print(f'{"sub":8s} {"contrast":11s} {"observed":>10s} {"null mean":>10s} '
          f'{"null sd":>9s} {"z":>7s}')
    allz = {c: [] for c in CONTRASTS}
    for s in subs:
        r = test(s, kpc)
        for c, (o, m, sd, z, n) in r.items():
            print(f'sub-{s:03d} {c:11s} {o:10.4f} {m:10.4f} {sd:9.4f} {z:7.2f}')
            allz[c].append(z)
        print(flush=True)
    print('median z across subjects: ' +
          '  '.join(f'{c} {np.median(v):+.2f}' for c, v in allz.items()))
    print(f'(|z| > 2 on a majority of subjects would be a real effect; '
          f'n={len(subs)} subjects)')
