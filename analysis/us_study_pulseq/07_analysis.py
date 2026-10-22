#!/usr/bin/env python
"""Statistics for the pulseq undersampling study (15 subjects x 4 gaze bins).

Inputs (results/us_pulseq/)
  metrics/region_metrics.csv          05_region_metrics.py (fixed masks)
  mriqc_{regu,grid}/group_T1w.tsv     02_run_mriqc.sh (MRIQC's own masks)
  metrics/mreyeqc_iqms_fixed.csv      06_mreyeqc_iqms.py iqms
  metrics/mreyeqc_ratings_astrained.csv  06_mreyeqc_iqms.py predict --labels as-trained
  readout_counts.csv                  readouts per subject x bin x level

Outputs (results/us_pulseq/analysis/)
  long.csv              every metric, one row per volume (x region)
  subject_level.csv     mean over the 4 gaze bins, plus the worst bin
  level_tests.csv       per method x region x metric x level: mean, SD, % change
                        vs clean, Wilcoxon (Holm within family), TOST +-5 %
  dose_response.csv     mixed models: metric ~ log2(readouts) vs ~ level

Unit of inference is the subject (n = 15): gaze bins are averaged first, so the
four bins of one subject never count as independent samples.
"""
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
OUT = ROOT / "analysis"
POS_IDX = {"up": 0, "down": 1, "left": 2, "right": 3}
LVL_MASK = {"full": "clean", "50p": "clean_0.50", "75p": "clean_0.75", "95p": "clean_0.95"}
LEVELS = ["woBin", "full", "50p", "75p", "95p"]
N_WOBIN = 79894
NYQUIST = np.pi * 240 ** 2 / 2             # full-diameter 3D radial spokes for N = 240

MRIQC_COLS = ["cjv", "cnr", "efc", "fber", "fwhm_avg", "inu_med", "qi_2", "snr_total", "snr_csf",
              "snr_gm", "snr_wm", "snrd_total", "snrd_csf", "snrd_gm", "snrd_wm", "wm2max",
              "icvs_csf", "icvs_gm", "icvs_wm", "tpm_overlap_csf", "tpm_overlap_gm", "tpm_overlap_wm"]
MREYEQC_COLS = ["seg_cnr", "seg_cjv", "seg_snr_LENS", "seg_snr_GLOBE", "seg_snr_OPTIQUE_NERVE",
                "seg_snr_MUSCLE", "seg_snr_FAT", "rpve_globe_boundary", "rpve_lens_boundary",
                "rpve_globe_lens_interface", "globe_to_max_intensity", "lens_to_max_intensity",
                "filter_laplace", "filter_sobel", "shannon_entropy", "ncc_window", "ssim_window",
                "nmi_window", "seg_globe_sphericity", "seg_lens_aspect_ratio"]

# +1: higher is better, -1: lower is better, 0: no direction
DIRECTION = {
    "ssim": 1, "ssim_wobin": 1, "nrmse": -1, "nrmse_wobin": -1, "dcte": 1, "dcte_crop": 1, "cnr": 1, "cjv": -1,
    "contrast": 0, "cv": -1, "edge_width": -1, "mean": 0, "std": 0, "p05": 0, "p95": 0,
    "efc": -1, "fber": 1, "fwhm_avg": -1, "inu_med": -1, "qi_2": -1, "wm2max": 1,
    "snr_total": 1, "snr_csf": 1, "snr_gm": 1, "snr_wm": 1,
    "snrd_total": 1, "snrd_csf": 1, "snrd_gm": 1, "snrd_wm": 1,
    "seg_cnr": 1, "seg_cjv": -1, "rating": 1,
}


def parse_bids(name):
    m = re.match(r"(sub-\d+)(?:_acq-(\w+?))?_rec-(steva|mathilda)(\w+?)_T1w", name)
    sub, pos, tag, lvl = m.groups()
    return dict(subject=sub, position=pos or "none", method={"steva": "REGU", "mathilda": "GRID"}[tag],
                level="woBin" if lvl == "Full" else lvl)


def readouts():
    rc = pd.read_csv(ROOT / "readout_counts.csv")
    rc["subject"] = rc["sub"].map(lambda s: f"sub-{s:03d}")
    inv_pos = {v: k for k, v in POS_IDX.items()}
    inv_lvl = {v: k for k, v in LVL_MASK.items()}
    rc["position"] = rc.region.map(inv_pos)
    rc["level"] = rc.lvl.map(inv_lvl)
    return rc[["subject", "position", "level", "used"]].rename(columns={"used": "n_readouts"})


def build_long():
    parts = []
    rm = pd.read_csv(ROOT / "metrics/region_metrics.csv")
    rm["source"] = "fixed"
    parts.append(rm)
    f = ROOT / "metrics/dcte_crop.csv"          # 05b_dcte_crop.py: DCTE on a box, no zero-filling
    if f.exists():
        dc = pd.read_csv(f)
        dc["source"] = "fixed"
        parts.append(dc)

    mq_all = [pd.read_csv(ROOT / f"mriqc_{t}/group_T1w.tsv", sep="\t")
              for t in ("regu", "grid", "wobin") if (ROOT / f"mriqc_{t}/group_T1w.tsv").exists()]
    if mq_all:
        # mriqc_wobin is a separate pass on the 15 woBin REGU volumes; it only fills
        # volumes whose IQM node crashed in the main run (sub-015 woBin did there)
        mq = pd.concat(mq_all, ignore_index=True).drop_duplicates("bids_name", keep="first")
        meta = pd.DataFrame([parse_bids(n) for n in mq.bids_name])
        mq = pd.concat([meta, mq[MRIQC_COLS].reset_index(drop=True)], axis=1).melt(
            id_vars=list(meta.columns), var_name="metric", value_name="value")
        mq["region"], mq["source"] = "mriqc", "mriqc"
        parts.append(mq)

    f = ROOT / "metrics/mreyeqc_iqms_fixed.csv"
    if f.exists():
        me = pd.read_csv(f)
        cols = [c for c in MREYEQC_COLS if c in me.columns]
        me = me.groupby(["subject", "method", "position", "level"])[cols].mean().reset_index()
        me = me.melt(id_vars=["subject", "method", "position", "level"], var_name="metric", value_name="value")
        me["region"], me["source"] = "mreyeqc", "mreyeqc"
        parts.append(me)

    f = ROOT / "metrics/mreyeqc_ratings_astrained.csv"
    if f.exists():
        ra = pd.read_csv(f).groupby(["subject", "method", "position", "level"]).mreyeqc_rating.mean()
        ra = ra.reset_index().rename(columns={"mreyeqc_rating": "value"})
        ra["metric"], ra["region"], ra["source"] = "rating", "mreyeqc", "mreyeqc_model"
        parts.append(ra)

    df = pd.concat(parts, ignore_index=True)
    df = df.merge(readouts(), on=["subject", "position", "level"], how="left")
    df.loc[df.level == "woBin", "n_readouts"] = N_WOBIN
    return df


def subject_level(df):
    """Mean over gaze bins, and the worst bin in the metric's own direction."""
    keys = ["source", "method", "region", "metric", "level", "subject"]
    sign = df.metric.map(DIRECTION).fillna(0).replace(0, np.nan)   # no direction -> no worst bin
    g = df.assign(_signed=df.value * sign).groupby(keys)
    out = pd.DataFrame({"mean_bins": g.value.mean(),
                        "worst_bin": g._signed.min() * g._signed.min().index.get_level_values("metric").map(DIRECTION),
                        "n_readouts": g.n_readouts.mean()})
    return out.reset_index()


def tost_paired(rel, margin=0.05):
    """Two one-sided t-tests on paired relative changes; equivalent if p < 0.05."""
    rel = np.asarray(rel, float)
    rel = rel[np.isfinite(rel)]
    if len(rel) < 3:
        return np.nan
    p_lo = stats.ttest_1samp(rel, -margin, alternative="greater").pvalue
    p_hi = stats.ttest_1samp(rel, margin, alternative="less").pvalue
    return max(p_lo, p_hi)


def level_tests(sl):
    rows = []
    for (src, method, region, metric), g in sl.groupby(["source", "method", "region", "metric"]):
        w = g.pivot(index="subject", columns="level", values="mean_bins")
        ref = "full" if "full" in w and w["full"].notna().sum() > 2 else None
        if ref is None and "woBin" in w:          # ssim / nrmse: reference is the clean bin itself
            ref = "woBin"
        fam = []
        for lvl in [l for l in LEVELS if l in w]:
            x = w[lvl].dropna()
            r = dict(source=src, method=method, region=region, metric=metric, level=lvl,
                     n=len(x), mean=x.mean(), sd=x.std(), median=x.median())
            if ref and lvl != ref:
                pair = w[[ref, lvl]].dropna()
                rel = (pair[lvl] - pair[ref]) / pair[ref].abs()
                r["pct_change"] = 100 * (pair[lvl].mean() - pair[ref].mean()) / abs(pair[ref].mean())
                r["pct_change_median_subj"] = 100 * rel.median()
                r["p_wilcoxon"] = stats.wilcoxon(pair[lvl], pair[ref]).pvalue if len(pair) > 5 else np.nan
                r["p_tost5"] = tost_paired(rel)
                fam.append(len(rows))
            rows.append(r)
        ps = [rows[i]["p_wilcoxon"] for i in fam]
        if ps and np.all(np.isfinite(ps)):
            adj = multipletests(ps, method="holm")[1]
            for i, p in zip(fam, adj):
                rows[i]["p_holm"] = p
    return pd.DataFrame(rows)


def dose_response(df, metrics):
    """Does quality follow the absolute readout count or the undersampling fraction?

    On every subject x bin x level volume (woBin excluded):
      A  value ~ log2(readouts) + (1 | subject)
      B  value ~ C(level)       + (1 | subject)
      C  value ~ log2(readouts) + C(level) + (1 | subject)
    Fitted by ML so AICs compare. If A beats B and C adds little, quality is set
    by how many readouts a bin has, whatever fraction of its clean data that is."""
    rows = []
    for (src, method, region, metric) in metrics:
        d = df[(df.source == src) & (df.method == method) & (df.region == region) & (df.metric == metric)
               & (df.level != "woBin")].dropna(subset=["value", "n_readouts"]).copy()
        if metric.startswith("ssim") and not metric.endswith("wobin"):
            d = d[d.level != "full"]                      # ssim(full, full) = 1 by construction
        d["l2n"] = np.log2(d.n_readouts)
        fits = {}
        for k, f in {"A": "value ~ l2n", "B": "value ~ C(level)", "C": "value ~ l2n + C(level)"}.items():
            try:
                fits[k] = smf.mixedlm(f, d, groups=d["subject"]).fit(reml=False)
            except Exception:
                fits[k] = None
        if fits["A"] is None:
            continue
        A = fits["A"]
        lr = 2 * (fits["C"].llf - A.llf) if fits["C"] is not None else np.nan
        dfree = d.level.nunique() - 1
        rows.append(dict(source=src, method=method, region=region, metric=metric, n=len(d),
                         slope_per_doubling=A.params["l2n"], slope_ci_lo=A.conf_int().loc["l2n", 0],
                         slope_ci_hi=A.conf_int().loc["l2n", 1], p_slope=A.pvalues["l2n"],
                         aic_readouts=A.aic, aic_level=fits["B"].aic if fits["B"] is not None else np.nan,
                         lr_level_given_readouts=lr, p_level_given_readouts=stats.chi2.sf(lr, dfree)))
    return pd.DataFrame(rows)


def fit_knee(l2n, y, subj, grid):
    """Profile least squares over the knee position k (subject fixed effects).

    value = a_subject + b * min(log2 N - k, 0): flat above k readouts, linear in
    log2 N below. Returns k, slope b, and the residual SS of the best and of the
    plain linear fit (for an F-type comparison). Subject intercepts are profiled
    out by within-subject centring, which is exact for fixed effects."""
    codes = pd.factorize(subj)[0]
    cnt = np.bincount(codes)

    def centre(v):
        return v - (np.bincount(codes, weights=v) / cnt)[codes]

    yc = centre(y)
    best = (np.inf, None, None)
    for k in grid:
        h = centre(np.minimum(l2n - k, 0.0))
        hh = h @ h
        b = (h @ yc) / hh if hh > 0 else 0.0
        rss = float(((yc - b * h) ** 2).sum())
        if rss < best[0]:
            best = (rss, k, b)
    lc = centre(l2n)
    bl = (lc @ yc) / (lc @ lc)
    rss_lin = float(((yc - bl * lc) ** 2).sum())
    return best[1], best[2], best[0], rss_lin


def knee_analysis(df, metrics, n_boot=2000, seed=0):
    """Readout count above which a metric stops improving (REGU saturates).

    Fitted on every subject x bin x level volume (woBin excluded). The knee's
    95 % CI resamples subjects with replacement. A knee at the edge of the grid
    means no plateau inside the sampled range."""
    rng = np.random.default_rng(seed)
    grid = np.arange(8.0, 14.01, 0.05)
    rows = []
    for (src, method, region, metric) in metrics:
        d = df[(df.source == src) & (df.method == method) & (df.region == region) & (df.metric == metric)
               & (df.level != "woBin")].dropna(subset=["value", "n_readouts"])
        if metric.startswith(("ssim", "nrmse")) and not metric.endswith("wobin"):
            continue                                    # relative to the bin's own clean volume
        l2n, y, subj = np.log2(d.n_readouts.values), d.value.values, d.subject.values
        if len(y) < 20 or np.nanstd(y) == 0:
            continue
        k, b, rss, rss_lin = fit_knee(l2n, y, subj, grid)
        n, p_lin = len(y), pd.get_dummies(subj).shape[1] + 1
        F = ((rss_lin - rss) / 1) / (rss / (n - p_lin - 1))
        subs = np.unique(subj)
        boots = []
        for _ in range(n_boot):
            pick = rng.choice(subs, len(subs), replace=True)
            idx = np.concatenate([np.flatnonzero(subj == s_) for s_ in pick])
            lab = np.concatenate([[f"{s_}_{i}"] * (subj == s_).sum() for i, s_ in enumerate(pick)])
            boots.append(fit_knee(l2n[idx], y[idx], lab, grid)[0])
        lo, hi = np.percentile(boots, [2.5, 97.5])
        plateau = float(np.mean(y[l2n >= k])) if (l2n >= k).any() else np.nan
        rows.append(dict(source=src, method=method, region=region, metric=metric, n=n,
                         knee_log2=k, knee_readouts=2 ** k, knee_ci_lo=2 ** lo, knee_ci_hi=2 ** hi,
                         knee_seconds=2 ** k * 0.00801, slope_below_per_doubling=b,
                         plateau_mean=plateau, F_vs_linear=F, p_vs_linear=stats.f.sf(F, 1, n - p_lin - 1),
                         at_grid_edge=bool(k <= grid[0] + 1e-9 or k >= grid[-1] - 1e-9)))
        print(f"  knee {method} {region} {metric}: {2 ** k:.0f} [{2 ** lo:.0f}, {2 ** hi:.0f}]", flush=True)
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(exist_ok=True)
    df = build_long()
    df.to_csv(OUT / "long.csv", index=False)
    sl = subject_level(df)
    sl.to_csv(OUT / "subject_level.csv", index=False)
    lt = level_tests(sl)
    lt.to_csv(OUT / "level_tests.csv", index=False)

    eye = [("fixed", m, r, k) for m in ("REGU", "GRID")
           for r, k in [("eyeball", "ssim"), ("lens", "ssim"), ("nerve", "ssim"), ("muscle", "ssim"),
                        ("eyes_mni", "ssim"), ("brain", "ssim"), ("wm", "ssim"), ("gm", "ssim"),
                        ("eyeball", "ssim_wobin"), ("brain", "ssim_wobin"), ("wm", "ssim_wobin"),
                        ("eyeball", "edge_width"), ("lens", "edge_width"),
                        ("lens_vitreous", "cnr"), ("nerve_fat", "cnr"), ("muscle_fat", "cnr"),
                        ("vitreous", "cv"), ("wm", "cv"), ("eyeball", "dcte"), ("brain", "dcte"),
                        ("eyeball", "dcte_crop"), ("brain", "dcte_crop")]]
    eye += [("mriqc", m, "mriqc", k) for m in ("REGU", "GRID") for k in ("cjv", "cnr", "efc", "fwhm_avg", "snr_wm")]
    eye += [("mreyeqc", m, "mreyeqc", k) for m in ("REGU", "GRID")
            for k in ("seg_cnr", "seg_cjv", "seg_snr_GLOBE", "seg_snr_MUSCLE", "seg_snr_OPTIQUE_NERVE", "filter_laplace")]
    eye += [("mreyeqc_model", m, "mreyeqc", "rating") for m in ("REGU", "GRID")]
    dr = dose_response(df, eye)
    dr.to_csv(OUT / "dose_response.csv", index=False)
    kn = knee_analysis(df, eye)
    kn.to_csv(OUT / "knee.csv", index=False)
    print(f"long {len(df)}  subject-level {len(sl)}  tests {len(lt)}  dose-response {len(dr)}")


if __name__ == "__main__":
    main()
