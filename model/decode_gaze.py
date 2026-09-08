"""Definitive gaze-decoding test: window sweep, both bin types, multi-shift null.

Supersedes probe_gaze.py / probe_gaze2.py, which were uninterpretable for two
reasons now fixed:

  * they scored plain accuracy, whose chance level is the retained-window
    majority rate (~0.33 here), not 0.25. Balanced accuracy is prior-invariant,
    so chance is 0.25 by construction.
  * they compared against a single circular shift. A *distribution* of shifts
    gives a null mean and spread, so the result is a z-score rather than a
    number to eyeball.

Everything else is as before: per-direction-cluster templates (the spoke
direction must be conditioned on), per-cluster local detrend, 32 PCs,
time-blocked CV, unsupervised centring of the test fold.
"""
import sys
import time
import numpy as np
import scipy.io as sio
import glob
from pathlib import Path
from scipy.spatial import cKDTree
from dprime import prep, CIDX, OUT, NR, NC

NFOLD = 5
NSEG, NSHOT, NOFF = 44, 1872, 14
WINDOWS = (64, 128, 256, 384, 640, 1024)      # 0.5 s .. 8 s at TR 8 ms
SHIFTS = (6000, 12000, 19000, 26000, 33000, 40000, 47000, 54000)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'


def gaze_bins(sub, kind):
    out = []
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/{kind}/eMask_th*region{r}.mat')
        m = np.asarray(sio.loadmat(g[0])['eMaskN']).squeeze().reshape(NSHOT, NSEG)
        out.append(m[NOFF:, 1:].ravel().astype(bool))
    return np.stack(out)


def readout_scores(Z, g):
    fold = np.arange(NR) * NFOLD // NR
    S = np.zeros((NR, 4), np.float32)
    seen = np.zeros(NR, bool)
    for f in range(NFOLD):
        te = fold == f
        for idx in CIDX:
            tr_i, te_i = idx[~te[idx]], idx[te[idx]]
            if te_i.size == 0 or tr_i.size < 40:
                continue
            mu = Z[tr_i].mean(0)
            T, keep = [None] * 4, []
            for k in range(4):
                sel = tr_i[g[k][tr_i]]
                if sel.size >= 5:
                    T[k] = Z[sel].mean(0) - mu
                    keep.append(k)
            if len(keep) < 2:
                continue
            Zt = Z[te_i] - Z[te_i].mean(0)
            sc = np.stack([(Zt @ np.conj(T[k])).real - 0.5 * np.vdot(T[k], T[k]).real
                           for k in keep], 1)
            sc = (sc - sc.mean()) / (sc.std() + 1e-9)
            S[np.ix_(te_i, np.array(keep))] = sc
            seen[te_i] = True
    return S, seen


def balanced_acc(S, seen, g, W, stride=64, purity=0.7):
    lab = np.full(NR, -1, np.int8)
    for k in range(4):
        lab[g[k]] = k
    oh = np.stack([(lab == k) for k in range(4)], 1).astype(np.int32)
    cs = np.concatenate([np.zeros((1, 4)), np.cumsum(np.where(seen[:, None], S, 0), 0)])
    cl = np.concatenate([np.zeros((1, 4), np.int64), np.cumsum(oh, 0)])
    cn = np.concatenate([[0], np.cumsum(seen)])
    hit = np.zeros(4); tot = np.zeros(4)
    for a in range(0, NR - W + 1, stride):
        b = a + W
        if cn[b] - cn[a] < 0.5 * W:
            continue
        c = cl[b] - cl[a]
        n = c.sum()
        if n < 0.3 * W or c.max() < purity * n:
            continue
        yt = int(np.argmax(c)); yh = int(np.argmax(cs[b] - cs[a]))
        tot[yt] += 1; hit[yt] += (yt == yh)
    ok = tot > 0
    return float(np.mean(hit[ok] / tot[ok])) if ok.sum() >= 3 else np.nan, int(tot.sum())


def run(sub, kind, shifts=SHIFTS):
    Z = prep(sub)
    g = gaze_bins(sub, kind)
    S, seen = readout_scores(Z, g)
    res = {}
    for W in WINDOWS:
        res[W] = [balanced_acc(S, seen, g, W)[0]]
    for sh in shifts:
        gs = np.roll(g, sh, axis=1)
        Ss, sn = readout_scores(Z, gs)
        for W in WINDOWS:
            res[W].append(balanced_acc(Ss, sn, gs, W)[0])
    return res


if __name__ == '__main__':
    kinds = ['clean', 'filtered']
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    for kind in kinds:
        print(f'\n=== bins: {kind} === balanced 4-way accuracy, chance 0.25')
        print(f'{"sub":8s} ' + ' '.join(f'{W*0.008:>5.1f}s' for W in WINDOWS))
        Z = {}
        for s in subs:
            t0 = time.time()
            r = run(s, kind)
            true = {W: r[W][0] for W in WINDOWS}
            z = {}
            for W in WINDOWS:
                nul = np.array(r[W][1:], float)
                z[W] = (true[W] - np.nanmean(nul)) / (np.nanstd(nul, ddof=1) + 1e-9)
            print(f'sub-{s:03d} ' + ' '.join(f'{true[W]:6.3f}' for W in WINDOWS) +
                  '   z: ' + ' '.join(f'{z[W]:+5.1f}' for W in WINDOWS) +
                  f'   [{time.time()-t0:.0f}s]', flush=True)
            Z[s] = z
        med = {W: np.median([Z[s][W] for s in subs]) for W in WINDOWS}
        print('median z   ' + ' '.join(f'{med[W]:+5.1f}' for W in WINDOWS))
