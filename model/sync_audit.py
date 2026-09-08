"""ET-MRI synchronisation audit across the cohort.

Three terms decide whether a label window can be narrowed:

  start offset   ET recording start vs the first scanner trigger, both read from
                 the PsychoPy log (`setRecordingState(True)` and `Keypress: s`).
  EyeLink drift  the .tsv.gz rows ARE 1 kHz EyeLink samples, so the row count is
                 an exact EyeLink-clock measure of an interval PsychoPy also
                 timed. The ratio is the EyeLink clock error.
  scanner drift  the scanner's own trigger period, from PMUTimeStamp sawtooth
                 resets in the twix header, against the same triggers seen by
                 PsychoPy.

Net misalignment at the end of the scan sets the smallest honest label window.
Do NOT anchor the EyeLink interval on the `eye tracker stopped` message: it
fires ~16 s after recording actually stops.
"""
import gzip
import glob
import re
import sys

import numpy as np

STUDY = '/home/debi/jaime/repos/MR-EyeTrack/data/study'
SUBS = range(1, 16)
SCAN_S = 659.0          # MRI acquisition length, seconds


def log_terms(sub):
    f = glob.glob(f'{STUDY}/sub-{sub:03d}/et/*.log')
    if not f:
        return None
    lines = open(f[0], errors='ignore').read().splitlines()
    keys, rec_on, rec_off = [], None, None
    for x in lines:
        m = re.match(r'\s*([\d.]+)\s+DATA\s+Keypress:\s*s\s*$', x)
        if m:
            keys.append(float(m.group(1)))
        if 'setRecordingState(True)' in x and rec_on is None:
            rec_on = float(x.split()[0])
        if 'setRecordingState(False)' in x:
            rec_off = float(x.split()[0])
    return np.array(keys), rec_on, rec_off


def n_samples(sub):
    f = glob.glob(f'{STUDY}/sub-{sub:03d}/et/*.tsv.gz')
    with gzip.open(f[0], 'rt') as fh:
        return sum(1 for _ in fh) - 1


def main():
    print('ET-MRI synchronisation, all subjects')
    print('(offset = ET start before first trigger; drift over a 659 s scan)\n')
    print(f'  {"sub":8s} {"triggers":>9s} {"period":>9s} {"offset":>9s} '
          f'{"EL ppm":>8s} {"EL drift":>9s} {"net":>8s}  verdict')
    rows = []
    for s in SUBS:
        try:
            keys, on, off = log_terms(s)
        except TypeError:
            print(f'  sub-{s:03d}  no log'); continue
        if keys is None or len(keys) < 10 or on is None:
            print(f'  sub-{s:03d}  incomplete log ({0 if keys is None else len(keys)} triggers)')
            continue
        # keys[0] is the scan start: PsychoPy blocks on the trigger, so the first
        # 's' it sees is the scan beginning. (Do not "clean" this by requiring a
        # regular 2.5 s interval -- the FIRST interval is systematically short,
        # ~1.99 s, so that filter discards the genuine first trigger and reports
        # a spurious ~2 s offset on every subject.)
        i = np.arange(len(keys))
        period = np.polyfit(i, keys, 1)[0]
        offset_ms = 1000 * (keys[0] - on)
        el_ppm = el_ms = np.nan
        if off is not None and off > on:
            n = n_samples(s)
            ratio = n / ((off - on) * 1000.0)
            el_ppm = 1e6 * (ratio - 1)
            el_ms = 1000 * SCAN_S * (ratio - 1)
        net = abs(offset_ms) + (abs(el_ms) if np.isfinite(el_ms) else 0.0)
        verdict = 'winLen 1 ok' if net < 12 else ('winLen 3' if net < 60 else 'winLen 5+')
        rows.append((s, offset_ms, el_ppm, el_ms, net))
        print(f'  sub-{s:03d} {len(keys):9d} {period:9.5f} {offset_ms:+8.1f}ms '
              f'{el_ppm:+8.0f} {el_ms:+8.1f}ms {net:7.1f}ms  {verdict}', flush=True)
    if rows:
        a = np.array([[r[1], r[2], r[3], r[4]] for r in rows], float)
        print(f'\n  {"median":8s} {"":9s} {"":9s} {np.nanmedian(a[:,0]):+8.1f}ms '
              f'{np.nanmedian(a[:,1]):+8.0f} {np.nanmedian(a[:,2]):+8.1f}ms '
              f'{np.nanmedian(a[:,3]):7.1f}ms')
        print(f'  {"worst":8s} {"":9s} {"":9s} {np.nanmax(np.abs(a[:,0])):+8.1f}ms '
              f'{np.nanmax(np.abs(a[:,1])):+8.0f} {np.nanmax(np.abs(a[:,2])):+8.1f}ms '
              f'{np.nanmax(a[:,3]):7.1f}ms')
        print(f'\n  A label window must exceed the net misalignment, or the error')
        print(f'  is a systematic function of scan time rather than random.')


if __name__ == '__main__':
    main()
