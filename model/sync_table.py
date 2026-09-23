"""Per-subject ET-MRI timing corrections, for the binning to apply.

Produces, for every subject, two numbers:

    offset_ms   elapsed ET time at the instant of the first MRI readout
    ratio       EyeLink samples per millisecond of real time

so the binning should index the ET mask at

    et_index(k) = ratio * (offset_ms + TimeStamp_ms(k))

instead of `TimeStamp_ms(k)`, which assumes offset 0 and ratio 1.

The offset resolves an ambiguity the earlier audit could not. Taking the first
`Keypress: s` in the log fails when a stray pulse from a preceding sequence
precedes it (sub-007, 009, 014 each have exactly one), and aligning the trigger
trains cannot help because the train is periodic at 2.5 s -- any whole-period
shift fits equally well. `PMUTimeStamp` at the *first readout* breaks it: it is
the time since the preceding trigger, so it pins the scan start *within* one
period. Choosing the trigger that puts the scan start nearest the ET recording
start then gives a unique answer.
"""
import contextlib
import glob
import gzip
import io
import re
import sys

import numpy as np
import scipy.io as sio
import scipy.integrate as _si

_si.cumtrapz = getattr(_si, 'cumtrapz', _si.cumulative_trapezoid)
sys.path.insert(0, '/home/debi/jaime/repos/MR-EyeTrack/old_study/code/twixtools')
import twixtools  # noqa: E402

STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
OUT = '/home/debi/jaime/repos/MR-EyeTrack/data/study/sync_corrections.mat'


def from_log(sub):
    L = glob.glob(f'{STUDY}/sub-{sub:03d}/et/*.log')[0]
    keys, on, off = [], None, None
    for x in open(L, errors='ignore'):
        m = re.match(r'\s*([\d.]+)\s+DATA\s+Keypress:\s*s\s*$', x)
        if m:
            keys.append(float(m.group(1)))
        if 'setRecordingState(True)' in x and on is None:
            on = float(x.split()[0])
        if 'setRecordingState(False)' in x:
            off = float(x.split()[0])
    return np.array(keys), on, off


def from_twix(sub):
    f = glob.glob(f'{STUDY}/sub-{sub:03d}/rawdata/*_T1wLIBRE.dat')[0]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        tw = twixtools.read_twix(f)
    img = [m for m in tw[-1]['mdb'] if m.is_image_scan()]
    ts = np.array([int(m.mdh.TimeStamp) for m in img], np.int64)
    order = np.argsort(ts)
    pmu0 = int(img[order[0]].mdh.PMUTimeStamp)
    return pmu0 * 2.5, len(img)          # ms since the preceding trigger


def n_et_samples(sub):
    f = glob.glob(f'{STUDY}/sub-{sub:03d}/et/*.tsv.gz')[0]
    with gzip.open(f, 'rt') as fh:
        return sum(1 for _ in fh) - 1


def main():
    subs = [int(x) for x in sys.argv[1:]] or list(range(1, 16))
    rows = []
    print(f'  {"sub":8s} {"pmu0":>8s} {"scan start (PC)":>16s} {"offset":>10s} '
          f'{"ppm":>7s} {"stray?":>7s}')
    for s in subs:
        keys, on, off = from_log(s)
        pmu0_ms, nro = from_twix(s)
        # the scan began pmu0_ms after some trigger; pick the trigger that puts
        # that instant closest to the ET recording start (unique within 2.5 s)
        cand = keys + pmu0_ms / 1000.0
        j = int(np.argmin(np.abs(cand - on)))
        scan_start = cand[j]
        offset_ms = (scan_start - on) * 1000.0
        ratio = n_et_samples(s) / ((off - on) * 1000.0) if off and off > on else 1.0
        stray = 'yes' if j > 0 else 'no'
        rows.append((s, offset_ms, ratio))
        print(f'  sub-{s:03d} {pmu0_ms:8.1f} {scan_start:16.4f} {offset_ms:+9.1f}ms '
              f'{1e6*(ratio-1):+7.0f} {stray:>7s}', flush=True)
    sio.savemat(OUT, {
        'subjects': np.array([r[0] for r in rows], float),
        'offset_ms': np.array([r[1] for r in rows], float),
        'ratio': np.array([r[2] for r in rows], float)})
    print(f'\nwrote {OUT}')
    print('  et_index(k) = ratio * (offset_ms + TimeStamp_ms(k))')


if __name__ == '__main__':
    main()
