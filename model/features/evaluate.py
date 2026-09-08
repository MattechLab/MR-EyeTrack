"""Do the CBMS-2025 style features work on this cohort, and on which task?

Two tasks, deliberately separated, because they are not the same problem:

  motion   eye moving vs eye still  -- the paper's task. Motion corrupts a
           readout's phase directly, and that mechanism is independent of where
           the eye is pointing, so direction-blind summary statistics can see it.
  gaze     which of 4 static positions -- needs the eye's *position* resolved,
           a sub-millimetre structural difference. Summary statistics collapse
           the sample axis and are blind to the spoke direction, so there is a
           real reason to expect them to do worse here.

Protocol follows the paper: transitional readouts excluded, classes balanced by
undersampling, features standardised per subject, and evaluation
leave-subjects-out. Grouped splits matter more than usual: a blink spans 25-50
consecutive readouts, so a random split would put the same event in train and
test and inflate everything.
"""
import glob
import sys
import numpy as np
import scipy.io as sio
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score, recall_score, precision_score

FEAT = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/features')
STUDY = Path('/home/debi/jaime/repos/MR-EyeTrack/data/study')
NR, NSEG, NSHOT, NOFF = 79894, 44, 1872, 14


def _mask(sub, path):
    m = np.asarray(sio.loadmat(path)['eMaskN']).squeeze().reshape(NSHOT, NSEG)
    return m[NOFF:, 1:].ravel().astype(bool)


def labels(sub, task):
    if task == 'motion':
        b = STUDY / f'sub-{sub:03d}/recon/bins'
        fix = _mask(sub, b / 'fixation-ok/eMask_th0.75_winLen10.mat')
        mot = np.zeros(NR, bool)
        for n in ('saccade-only', 'blink-event', 'tracking-loss'):
            mot |= _mask(sub, b / n / 'eMask_th0.75_winLen10.mat')
        y = np.full(NR, -1, np.int8)
        y[fix & ~mot] = 0
        y[mot] = 1                       # readouts in neither are transitional: dropped
        return y
    y = np.full(NR, -1, np.int8)
    for r in range(4):
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/recon/bins/filtered/eMask_th*region{r}.mat')
        y[_mask(sub, g[0])] = r
    return y


def load(sub, task):
    X = np.load(FEAT / f'sub-{sub:03d}_readout_features.npy')
    y = labels(sub, task)
    keep = y >= 0
    X = X[keep]; y = y[keep]
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)
    X = StandardScaler().fit_transform(X)        # per subject, as the paper does
    return X.astype(np.float32), y, np.flatnonzero(keep)


def balance(y, rng, groups=None):
    """Undersample the majority class; keep whole events together if grouped."""
    idx = []
    n = min((y == c).sum() for c in np.unique(y))
    for c in np.unique(y):
        i = np.flatnonzero(y == c)
        idx.append(rng.choice(i, n, replace=False))
    return np.sort(np.concatenate(idx))


def main():
    task = sys.argv[1] if len(sys.argv) > 1 else 'motion'
    subs = [int(a) for a in sys.argv[2:]] or None
    avail = sorted(int(p.name[4:7]) for p in FEAT.glob('sub-*_readout_features.npy'))
    subs = subs or avail
    subs = [s for s in subs if s in avail]
    print(f'task={task}  subjects available: {subs}')
    if len(subs) < 3:
        print('need at least 3 subjects for leave-subjects-out'); return

    rng = np.random.default_rng(0)
    D = {s: load(s, task) for s in subs}
    for s in subs:
        _, y, _ = D[s]
        print(f'  sub-{s:03d}: ' + ', '.join(f'class {c} n={int((y==c).sum())}'
                                             for c in np.unique(y)))
    nclass = len(np.unique(np.concatenate([D[s][1] for s in subs])))
    chance = 1.0 / nclass
    print(f'\nleave-one-subject-out, balanced accuracy (chance {chance:.3f})')
    print(f'  {"held-out":10s} {"HistGB":>8s} {"LinearSVC":>10s}')
    res = {'HistGB': [], 'LinearSVC': []}
    for held in subs:
        tr = [s for s in subs if s != held]
        Xtr = np.concatenate([D[s][0] for s in tr])
        ytr = np.concatenate([D[s][1] for s in tr])
        sel = balance(ytr, rng)
        Xtr, ytr = Xtr[sel], ytr[sel]
        Xte, yte, _ = D[held]
        row = {}
        for nm, clf in (('HistGB', HistGradientBoostingClassifier(max_iter=200,
                                                                 random_state=0)),
                        ('LinearSVC', LinearSVC(C=0.01, dual='auto', max_iter=3000))):
            clf.fit(Xtr, ytr)
            b = balanced_accuracy_score(yte, clf.predict(Xte))
            row[nm] = b; res[nm].append(b)
        print(f'  sub-{held:03d}    {row["HistGB"]:8.4f} {row["LinearSVC"]:10.4f}', flush=True)
    print(f'  {"median":10s} {np.median(res["HistGB"]):8.4f} '
          f'{np.median(res["LinearSVC"]):10.4f}')


if __name__ == '__main__':
    main()
