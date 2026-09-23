"""How much of each ET class survives into the readout labels?

The ET runs at 1 kHz and the sequence at TR = 8 ms, so **1 readout ~ 8 ET
samples**. If labelling were lossless, a class covering N ET samples would cover
about N/8 readouts.

It is not lossless. `eyeGenerateBinningWin` marks a readout only when 75% of a
winLen=10 window (80 ms) agrees, which is a purity filter designed for
reconstruction binning. Events shorter than 60 ms cannot reach that threshold at
any position, so saccades (30-80 ms) are erased while blinks (100-400 ms) pass.

This table makes the loss explicit: expected readouts (ET samples / 8) against
the readouts actually labelled.
"""
import glob
import sys

import numpy as np
import scipy.io as sio

STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
NSEG, NSHOT, NOFF = 44, 1872, 14
NACQ = NSEG * NSHOT                    # 82368 acquired readouts
MS_PER_READOUT = 8.0

# class -> (ET mask filename fragment, readout mask path fragment)
CLASSES = [
    ('fixation',      'fixation_ok_mask',   'fixation-ok/eMask_th0.75_winLen{w}.mat'),
    ('up',            'clean_mask_0.3_0',   'clean/eMask_th0.75_region0.mat|clean/eMask_th0.75_winLen{w}_region0.mat'),
    ('down',          'clean_mask_0.3_1',   'clean/eMask_th0.75_region1.mat|clean/eMask_th0.75_winLen{w}_region1.mat'),
    ('left',          'clean_mask_0.3_2',   'clean/eMask_th0.75_region2.mat|clean/eMask_th0.75_winLen{w}_region2.mat'),
    ('right',         'clean_mask_0.3_3',   'clean/eMask_th0.75_region3.mat|clean/eMask_th0.75_winLen{w}_region3.mat'),
    ('saccade',       'saccade_only_mask',  'saccade-only/eMask_th0.75_winLen{w}.mat'),
    ('blink',         'blink_event_mask',   'blink-event/eMask_th0.75_winLen{w}.mat'),
    ('tracking-loss', 'tracking_loss_mask', 'tracking-loss/eMask_th0.75_winLen{w}.mat'),
]
WINLENS = ('10', '3', '3_sync')


def load_mask(path):
    d = sio.loadmat(path)
    k = [x for x in d if not x.startswith('__')][0]
    return np.asarray(d[k]).squeeze().astype(bool)


def counts(sub):
    out = {}
    for name, etfrag, rofrag in CLASSES:
        n_et = None
        g = glob.glob(f'{STUDY}/sub-{sub:03d}/eyemasks/*{etfrag}.mat')
        if g:
            n_et = int(load_mask(g[0]).sum())
        per_w = {}
        for w in WINLENS:
            n_ro = None
            for cand in rofrag.format(w=w).split('|'):
                # the winLen10 4-way masks predate the winLen tag in the filename
                if w != '10' and 'winLen' not in cand:
                    continue
                if w == '10' and 'winLen10_region' in cand:
                    continue
                try:
                    n_ro = int(load_mask(f'{STUDY}/sub-{sub:03d}/recon/bins/{cand}').sum())
                    break
                except Exception:
                    continue
            per_w[w] = n_ro
        out[name] = (n_et, per_w)
    return out


def main():
    subs = [int(x) for x in sys.argv[1:]] or list(range(1, 16))
    print(f'ET samples (1 ms) vs readouts (8 ms).  1 readout ~ {MS_PER_READOUT:.0f} ET samples.')
    print('"expected" = ET samples / 8 ; "kept" = expected that survived the '
          'winLen=10 th=0.75 window\n')
    agg = {c[0]: [0, 0, 0, 0] for c in CLASSES}
    for s in subs:
        c = counts(s)
        print(f'  sub-{s:03d}')
        print(f'    {"class":14s} {"ET samples":>11s} {"expect":>8s} '
              f'{"w10":>7s} {"w3":>7s} {"w3+sync":>7s}')
        for name, _, _ in CLASSES:
            n_et, per_w = c[name]
            if n_et is None:
                continue
            exp = n_et / MS_PER_READOUT
            a, b, c3 = per_w.get('10'), per_w.get('3'), per_w.get('3_sync')
            f = lambda v: f'{100*v/exp:6.1f}%' if v is not None and exp > 0 else '     -'
            agg[name][0] += n_et
            agg[name][1] += a or 0
            agg[name][2] += b or 0
            agg[name][3] += c3 or 0
            print(f'    {name:14s} {n_et:11d} {exp:8.0f} {f(a)} {f(b)} {f(c3)}')
        print()
    print('  COHORT TOTAL')
    print(f'    {"class":14s} {"expect":>9s} {"w10":>16s} {"w3":>16s} {"w3+sync":>16s}')
    for name, _, _ in CLASSES:
        n_et, a, b, c3 = agg[name]
        exp = n_et / MS_PER_READOUT
        g = lambda v: f'{v:8d} {100*v/exp:5.1f}%' if exp else '        -'
        print(f'    {name:14s} {exp:9.0f} {g(a):>16s} {g(b):>16s} {g(c3):>16s}')


if __name__ == '__main__':
    main()
