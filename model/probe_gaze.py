"""Linear positive control: is gaze state decodable from a window of readouts?

Nearest-mean classifier in the per-direction-cluster residual space, which is
the simplest thing that respects the geometry: a readout's profile depends on
its spoke direction, so class templates are learned per direction cluster and a
readout is only ever compared against templates for its own cluster.

Guards against the ways this can fool you:
  * time-blocked CV -- templates never see the fold they score
  * test folds are centred on their OWN cluster means (unsupervised), so a
    slow drift between train and test cannot masquerade as class evidence
  * scores standardised per cluster before pooling, so window predictions
    cannot ride on which clusters a window happens to contain
  * a circular-shift null re-runs everything with the labels rotated in time
"""
import sys
import numpy as np
from pathlib import Path
from scipy.spatial import cKDTree

OUT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/kband')
NC, NFOLD = 600, 5
WINDOWS = (32, 128, 256, 512, 1024)
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)

dirs = np.load(OUT / 'dirs.npy')
NR = dirs.shape[0]
rng = np.random.default_rng(0)
CL = cKDTree(dirs[rng.choice(NR, NC, replace=False)]).query(dirs)[1]
CIDX = [np.where(CL == c)[0] for c in range(NC)]


def readout_scores(X, gaze, shift=0):
    """Per-readout, per-class evidence, from time-blocked CV. [NR, 4]"""
    g = np.roll(gaze, shift, axis=1) if shift else gaze
    fold = np.arange(NR) * NFOLD // NR
    S = np.full((NR, 4), np.nan, np.float32)
    for f in range(NFOLD):
        te = fold == f
        for c in range(NC):
            idx = CIDX[c]
            tr_i, te_i = idx[~te[idx]], idx[te[idx]]
            if te_i.size == 0 or tr_i.size < 40:
                continue
            mu_tr = X[tr_i].mean(0)
            T, keep = [], []
            for k in range(4):
                sel = tr_i[g[k][tr_i]]
                if sel.size < 5:
                    T.append(None); continue
                T.append(X[sel].mean(0) - mu_tr); keep.append(k)
            if len(keep) < 2:
                continue
            Xt = X[te_i] - X[te_i].mean(0)          # unsupervised test centring
            for k in keep:
                S[te_i, k] = (Xt @ np.conj(T[k])).real - 0.5 * np.vdot(T[k], T[k]).real
            kk = np.array(keep)
            blk = S[te_i][:, kk]
            S[np.ix_(te_i, kk)] = (blk - blk.mean()) / (blk.std() + 1e-9)
    return S


def window_eval(S, gaze, W, stride=32, purity=0.6, min_lab=0.3):
    """Pool scores over W consecutive readouts and predict one label per window."""
    ok = ~np.isnan(S).any(1)
    Sf = np.where(np.isnan(S), 0.0, S)
    cs = np.concatenate([np.zeros((1, 4)), np.cumsum(Sf, 0)])
    cn = np.concatenate([[0], np.cumsum(ok)])
    lab = np.full(NR, -1, np.int8)
    for k in range(4):
        lab[gaze[k]] = k
    onehot = np.stack([(lab == k) for k in range(4)], 1).astype(np.int32)
    cl = np.concatenate([np.zeros((1, 4), np.int64), np.cumsum(onehot, 0)])
    yhat, ytrue = [], []
    for a in range(0, NR - W + 1, stride):
        b = a + W
        if cn[b] - cn[a] < 0.5 * W:
            continue
        counts = cl[b] - cl[a]
        n = counts.sum()
        if n < min_lab * W or counts.max() < purity * n:
            continue
        yhat.append(int(np.argmax(cs[b] - cs[a])))
        ytrue.append(int(np.argmax(counts)))
    return np.array(ytrue), np.array(yhat)


def run(sub, shift=0):
    X = np.load(OUT / f'sub-{sub:03d}_band160-320_ROI-PCA_8.npy').reshape(NR, -1)
    X = X / (np.abs(X).mean() + 1e-12)
    gaze = np.load(OUT / f'sub-{sub:03d}_labels.npz')['gaze']
    S = readout_scores(X, gaze, shift)
    out = {}
    for W in WINDOWS:
        yt, yh = window_eval(S, gaze, W)
        out[W] = (float((yt == yh).mean()) if yt.size else np.nan, int(yt.size))
    return out


if __name__ == '__main__':
    subs = [int(a) for a in sys.argv[1:]] or list(KEEP)
    print('4-way gaze accuracy from a window of readouts (chance 0.25), '
          f'within subject, {NFOLD}-fold time-blocked CV')
    print('window   ' + ' '.join(f'{W*0.008:5.1f}s' for W in WINDOWS) + '\n')
    acc = {W: [] for W in WINDOWS}
    for s in subs:
        r = run(s)
        print(f'sub-{s:03d} ' + ' '.join(f'{r[W][0]:6.3f}' for W in WINDOWS) +
              '   n=' + ' '.join(f'{r[W][1]}' for W in WINDOWS), flush=True)
        for W in WINDOWS:
            acc[W].append(r[W][0])
    print('median  ' + ' '.join(f'{np.nanmedian(acc[W]):6.3f}' for W in WINDOWS))
    print('\ncircular-shift null (labels rotated 20 000 readouts = 160 s):')
    nul = {W: [] for W in WINDOWS}
    for s in subs[:4]:
        r = run(s, shift=20000)
        print(f'sub-{s:03d} ' + ' '.join(f'{r[W][0]:6.3f}' for W in WINDOWS), flush=True)
        for W in WINDOWS:
            nul[W].append(r[W][0])
    print('median  ' + ' '.join(f'{np.nanmedian(nul[W]):6.3f}' for W in WINDOWS))
