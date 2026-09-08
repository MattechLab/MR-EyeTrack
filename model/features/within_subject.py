"""Control: is the motion information absent, or just subject-specific?

Leave-one-subject-out came out at exactly chance. Two very different causes:

  (a) the compressed data carries no motion information at all, or
  (b) it carries it, but in a per-subject basis that cannot transfer.

(b) is a live hypothesis because ROI-PCA is estimated per subject -- virtual
channel 3 of sub-001 has no relationship to sub-004's. The CBMS paper used 20
*physical* coils, which mean the same thing in every subject.

A within-subject time-blocked split separates the two. Time-blocked, not random:
a blink spans 25-50 consecutive readouts, so a random split puts the same event
on both sides and inflates everything.
"""
import numpy as np
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import balanced_accuracy_score
from evaluate import labels, FEAT, NR

KEEP = (1, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15)
rng = np.random.default_rng(0)
print('motion vs still, WITHIN subject, time-blocked 70/30 (chance 0.500)')
print(f'  {"sub":8s} {"train n":>9s} {"test n":>8s} {"balanced acc":>13s}')
out = []
for s in KEEP:
    X = np.nan_to_num(np.load(FEAT / f'sub-{s:03d}_readout_features.npy'),
                      nan=0., posinf=0., neginf=0.)
    y = labels(s, 'motion')
    cut = int(0.7 * NR)
    tr = (np.arange(NR) < cut) & (y >= 0)
    te = (np.arange(NR) >= cut) & (y >= 0)
    if (y[tr] == 1).sum() < 50 or (y[te] == 1).sum() < 50:
        print(f'  sub-{s:03d}  too few motion readouts'); continue
    sc = StandardScaler().fit(X[tr])
    Xtr, ytr = sc.transform(X[tr]), y[tr]
    Xte, yte = sc.transform(X[te]), y[te]
    n = min((ytr == 0).sum(), (ytr == 1).sum())
    sel = np.concatenate([rng.choice(np.flatnonzero(ytr == c), n, replace=False)
                          for c in (0, 1)])
    clf = HistGradientBoostingClassifier(max_iter=200, random_state=0)
    clf.fit(Xtr[sel], ytr[sel])
    b = balanced_accuracy_score(yte, clf.predict(Xte))
    out.append(b)
    print(f'  sub-{s:03d} {2*n:9d} {int(te.sum()):8d} {b:13.4f}', flush=True)
print(f'  {"median":8s} {"":9s} {"":8s} {np.median(out):13.4f}')
print('\n  at chance within subject too -> the information is not there.')
print('  well above chance within subject -> it is there but subject-specific,')
print('  which is what a per-subject ROI-PCA basis would cause.')
