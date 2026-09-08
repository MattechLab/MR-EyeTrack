"""Same tests as evaluate.py, on the uncompressed 52-channel features.

Three tests, in the order that makes the result interpretable:

  1. within-subject motion (time-blocked)  -- is the information there at all?
  2. leave-one-subject-out motion          -- the CBMS protocol; does it transfer?
  3. within-subject gaze (4-class)         -- does uncompressed help gaze too?

Test 1 first, because on the compressed data both 1 and 2 were at chance, which
proved the information was absent rather than merely untransferable. If 1 is
above chance here and the compressed version was not, the compression is the
cause, and that is the whole point of this run.

Rows are pulled from a memory-mapped file so only the selected readouts are
read -- the full matrix is 1.5 GiB per subject.
"""
import sys
import numpy as np
from pathlib import Path
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score
from evaluate import labels as _labels
from extract_raw import NR

RAW = Path('/home/debi/jaime/repos/MR-EyeTrack/data/derived/features_raw')
KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)


def available():
    return sorted(int(p.name[4:7]) for p in RAW.glob('sub-*_raw_features.npy'))


def rows(sub, idx):
    X = np.load(RAW / f'sub-{sub:03d}_raw_features.npy', mmap_mode='r')
    return np.nan_to_num(np.asarray(X[idx]), nan=0., posinf=0., neginf=0.)


def gaze_labels(sub):
    return _labels(sub, 'gaze')


def motion_labels(sub):
    return _labels(sub, 'motion')


def balanced_fit_predict(Xtr, ytr, Xte, yte, rng, seed=0):
    n = min((ytr == c).sum() for c in np.unique(ytr))
    sel = np.concatenate([rng.choice(np.flatnonzero(ytr == c), n, replace=False)
                          for c in np.unique(ytr)])
    sc = StandardScaler().fit(Xtr[sel])
    clf = HistGradientBoostingClassifier(max_iter=200, random_state=seed)
    clf.fit(sc.transform(Xtr[sel]), ytr[sel])
    return balanced_accuracy_score(yte, clf.predict(sc.transform(Xte))), 2 * n


def within(subs, task):
    rng = np.random.default_rng(0)
    lab = motion_labels if task == 'motion' else gaze_labels
    nclass = 2 if task == 'motion' else 4
    print(f'\n=== within-subject, time-blocked 70/30, {task} '
          f'(chance {1/nclass:.3f}) ===')
    print(f'  {"sub":8s} {"train n":>9s} {"test n":>8s} {"balanced acc":>13s}')
    out = []
    for s in subs:
        y = lab(s)
        cut = int(0.7 * NR)
        tr = np.flatnonzero((np.arange(NR) < cut) & (y >= 0))
        te = np.flatnonzero((np.arange(NR) >= cut) & (y >= 0))
        if min(np.bincount(y[tr], minlength=nclass)) < 50:
            print(f'  sub-{s:03d}  too few in one class'); continue
        b, ntr = balanced_fit_predict(rows(s, tr), y[tr], rows(s, te), y[te], rng)
        out.append(b)
        print(f'  sub-{s:03d} {ntr:9d} {len(te):8d} {b:13.4f}', flush=True)
    if out:
        print(f'  {"median":8s} {"":9s} {"":8s} {np.median(out):13.4f}')
    return out


def loso(subs, task):
    rng = np.random.default_rng(0)
    lab = motion_labels if task == 'motion' else gaze_labels
    nclass = 2 if task == 'motion' else 4
    print(f'\n=== leave-one-subject-out, {task} (chance {1/nclass:.3f}) ===')
    print(f'  {"held-out":10s} {"balanced acc":>13s}')
    cache = {}
    for s in subs:
        y = lab(s)
        k = np.flatnonzero(y >= 0)
        cache[s] = (k, y[k])
    out = []
    for held in subs:
        tr = [s for s in subs if s != held]
        # subsample per training subject to keep memory sane
        Xs, ys = [], []
        for s in tr:
            k, yy = cache[s]
            n = min(len(k), 20000)
            pick = rng.choice(len(k), n, replace=False)
            Xs.append(rows(s, k[pick])); ys.append(yy[pick])
        Xtr = np.concatenate(Xs); ytr = np.concatenate(ys)
        del Xs
        k, yte = cache[held]
        b, _ = balanced_fit_predict(Xtr, ytr, rows(held, k), yte, rng)
        out.append(b)
        print(f'  sub-{held:03d}    {b:13.4f}', flush=True)
    if out:
        print(f'  {"median":10s} {np.median(out):13.4f}')
    return out


if __name__ == '__main__':
    subs = [int(a) for a in sys.argv[1:]] or available()
    subs = [s for s in subs if s in available()]
    print(f'uncompressed 52-channel features, subjects: {subs}')
    print(f'compressed reference: motion within 0.5018, LOSO 0.4984 (both chance)')
    within(subs, 'motion')
    if len(subs) >= 3:
        loso(subs, 'motion')
    within(subs, 'gaze')
