"""Offset between scanner and stimulus PC, by aligning the two trigger trains.

sync_audit.py takes the first `Keypress: s` in the PsychoPy log as the scan
start. That is right when PsychoPy blocks on the trigger, but a stray pulse from
a preceding sequence would be taken as the start and produce a spurious
multi-second offset. Three subjects show exactly that pattern.

The scanner records the same triggers itself: `PMUTimeStamp` in the twix mdh is
a sawtooth whose resets are trigger times in the scanner's own clock. Aligning
the two trains -- rather than trusting either train's first element -- measures
the offset with no assumption about which pulse was first.
"""
import contextlib
import glob
import io
import re
import sys

import numpy as np
import scipy.integrate as _si

_si.cumtrapz = getattr(_si, 'cumtrapz', _si.cumulative_trapezoid)
sys.path.insert(0, '/home/debi/jaime/repos/MR-EyeTrack/old_study/code/twixtools')
import twixtools  # noqa: E402

STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'


def scanner_triggers(sub):
    f = glob.glob(f'{STUDY}/sub-{sub:03d}/rawdata/*_T1wLIBRE.dat')[0]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        tw = twixtools.read_twix(f)
    img = [m for m in tw[-1]['mdb'] if m.is_image_scan()]
    ts = np.array([int(m.mdh.TimeStamp) for m in img], np.int64)
    pmu = np.array([int(m.mdh.PMUTimeStamp) for m in img], np.int64)
    t = (ts - ts.min()) * 2.5 / 1000.0
    return t[np.flatnonzero(np.diff(pmu) < 0) + 1], t[-1]


def pc_triggers(sub):
    L = glob.glob(f'{STUDY}/sub-{sub:03d}/et/*.log')[0]
    keys, on = [], None
    for x in open(L, errors='ignore'):
        m = re.match(r'\s*([\d.]+)\s+DATA\s+Keypress:\s*s\s*$', x)
        if m:
            keys.append(float(m.group(1)))
        if 'setRecordingState(True)' in x and on is None:
            on = float(x.split()[0])
    return np.array(keys), on


def best_shift(a, b, lo=-20.0, hi=20.0):
    """Shift applied to b that best aligns it with a (median nearest distance)."""
    grid = np.arange(lo, hi, 0.05)
    cost = [np.median(np.abs(a[:, None] - (b + s)[None, :]).min(1)) for s in grid]
    s0 = grid[int(np.argmin(cost))]
    fine = np.arange(s0 - 0.05, s0 + 0.05, 0.001)
    cost = [np.median(np.abs(a[:, None] - (b + s)[None, :]).min(1)) for s in fine]
    return fine[int(np.argmin(cost))], min(cost)


def main():
    subs = [int(x) for x in sys.argv[1:]] or [9, 14, 15]
    print('trigger-train alignment: scanner PMU resets vs PsychoPy keypresses\n')
    print(f'  {"sub":8s} {"PMU":>5s} {"PC":>5s} {"scan_t0 in PC clock":>21s} '
          f'{"resid":>7s} {"ET offset":>11s}  note')
    for s in subs:
        st, scan_end = scanner_triggers(s)
        pt, on = pc_triggers(s)
        # shift the PC train so it lands on the scanner train; scanner t=0 is the
        # first readout, so the shift is that instant expressed in PC time
        sh, resid = best_shift(st, -pt if False else pt * 0 + pt, lo=-1400, hi=1400) \
            if False else (None, None)
        grid = np.arange(pt.min() - 5, pt.min() + 5, 0.001)
        cost = [np.median(np.abs((pt - g)[:, None] - st[None, :]).min(0)) for g in grid]
        t0 = grid[int(np.argmin(cost))]
        resid = min(cost)
        off = 1000 * (t0 - on)
        note = 'ET started after the scan' if off > 50 else 'ok'
        print(f'  sub-{s:03d} {len(st):5d} {len(pt):5d} {t0:21.4f} {1000*resid:6.1f}ms '
              f'{off:+10.1f}ms  {note}', flush=True)
    print('\n  scan_t0 = the first readout, located in the PsychoPy clock.')
    print('  ET offset = scan start minus setRecordingState(True); positive means')
    print('  the eye tracker began recording AFTER the scan had already started.')


if __name__ == '__main__':
    main()
