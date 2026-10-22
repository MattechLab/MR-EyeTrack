#!/usr/bin/env python
"""Figures for the pulseq undersampling study (abstract + supplementary).

    python 08_figures.py [fig1 fig2 ...]      # default: all

Writes PNG (300 dpi) and PDF to results/us_pulseq/figs/.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy import ndimage as ndi

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "results/us_pulseq"
BIDS = ROOT / "bids/derivatives/reconstructions"
AN = ROOT / "analysis"
FIGS = ROOT / "figs"

# reference palette (dataviz skill): categorical slots 1-2 for the two methods,
# one-hue blue ordinal ramp (steps 700 -> 250) for clean -> 95 % undersampling
C_REGU, C_GRID = "#2a78d6", "#eb6834"
LEVEL_RAMP = {"full": "#0d366b", "50p": "#1c5cab", "75p": "#3987e5", "95p": "#86b6ef"}
INK, INK2, GRIDC = "#0b0b0b", "#52514e", "#e4e3df"
LEVELS = ["woBin", "full", "50p", "75p", "95p"]
LABEL = {"woBin": "No binning", "full": "Clean", "50p": "50 %", "75p": "75 %", "95p": "95 %"}
TR_S = 0.00801                      # 659.77 s / 82 368 readouts
NYQ = np.pi * 240 ** 2 / 2

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRIDC,
    "grid.linewidth": 0.6, "legend.frameon": False, "lines.linewidth": 2, "figure.dpi": 100,
})


def save(fig, name):
    FIGS.mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("wrote", FIGS / f"{name}.png")


def vol(sub, method, pos, level):
    tag = {"REGU": "steva", "GRID": "mathilda"}[method]
    p = (BIDS / sub / "anat" / f"{sub}_rec-{tag}Full_T1w.nii.gz" if level == "woBin"
         else BIDS / sub / "anat" / f"{sub}_acq-{pos}_rec-{tag}{level}_T1w.nii.gz")
    return nib.load(p).get_fdata(dtype=np.float32)


def orbit_crop(sub):
    L = np.asanyarray(nib.load(ROOT / "masks" / sub / "aeye_labels.nii.gz").dataobj)
    z = int(round(ndi.center_of_mass(L == 1)[2]))
    xs = np.nonzero(L.any(axis=(1, 2)))[0]
    ys = np.nonzero((L == 2).any(axis=(0, 2)))[0]
    return L, z, (slice(xs.min() - 6, xs.max() + 7), slice(ys.min() - 26, ys.max() + 7))


def fig1(sub="sub-004", pos="left"):
    """A: example orbit + sagittal images per method and level. B: readouts per bin."""
    rc = pd.read_csv(ROOT / "readout_counts.csv")
    lv = {"clean": "full", "clean_0.50": "50p", "clean_0.75": "75p", "clean_0.95": "95p"}
    rc["level"] = rc.lvl.map(lv)
    L, z, (sx, sy) = orbit_crop(sub)
    pidx = {"up": 0, "down": 1, "left": 2, "right": 3}[pos]
    n_here = rc[(rc["sub"] == int(sub[-3:])) & (rc.region == pidx)].set_index("level").used
    orb_aspect = (sy.stop - sy.start) / (sx.stop - sx.start)
    sag_y, sag_z = slice(25, 225), slice(20, 220)

    fig = plt.figure(figsize=(7.2, 4.6))
    gsA = fig.add_gridspec(4, 5, left=0.04, right=0.74, top=0.88, bottom=0.02, wspace=0.03, hspace=0.04,
                           height_ratios=[orb_aspect, 1, orb_aspect, 1])
    for i, method in enumerate(["REGU", "GRID"]):
        for j, lvl in enumerate(LEVELS):
            v = vol(sub, method, pos, lvl)
            orb = v[sx, sy, z]
            sag = v[v.shape[0] // 2, sag_y, sag_z]
            for k, (img, vmax) in enumerate([(orb, np.percentile(orb, 99.5)), (sag, 1100)]):
                ax = fig.add_subplot(gsA[2 * i + k, j])
                ax.imshow(img.T, cmap="gray", origin="lower", vmin=0, vmax=vmax, interpolation="nearest",
                          aspect="auto")
                ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
                for sp in ax.spines.values():
                    sp.set_visible(False)
                if k == 0 and lvl == "woBin" and method == "REGU":
                    lab = L[sx, sy, z]
                    for val, c in [(1, "#eda100"), (2, "#1baf7a"), (3, "#e87ba4"), (6, "#e34948"), (7, "#e34948")]:
                        if (lab == val).any():
                            ax.contour((lab == val).T, levels=[0.5], colors=c, linewidths=0.6)
                if i == 0 and k == 0:
                    n = 79894 if lvl == "woBin" else int(n_here[lvl])
                    ax.set_title(f"{LABEL[lvl]}\n{n:,}", fontsize=7.5, color=INK, pad=3)
                if j == 0:
                    ax.set_ylabel(method if k == 0 else "", fontsize=8, color=INK, labelpad=2)
    fig.text(0.005, 0.975, "A", fontsize=10, weight="bold", color=INK)
    fig.text(0.39, 0.985, f"one subject, '{pos}' gaze bin: readouts used", fontsize=7, color=INK2, ha="center")

    ax = fig.add_axes([0.82, 0.12, 0.17, 0.78])
    rng = np.random.default_rng(0)
    order = ["full", "50p", "75p", "95p"]
    for x, lvl in enumerate(order):
        y = rc[rc.level == lvl].used.values
        ax.scatter(x + rng.uniform(-0.2, 0.2, len(y)), y, s=8, color=LEVEL_RAMP[lvl], alpha=0.85,
                   edgecolor="white", linewidth=0.3, zorder=3)
        ax.plot([x - 0.28, x + 0.28], [np.median(y)] * 2, color=INK, lw=1.5, zorder=4)
    ax.set_yscale("log", base=2)
    ticks = [128, 512, 2048, 8192, 32768]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t // 1024}k" if t >= 1024 else str(t) for t in ticks])
    ax.minorticks_off()
    for t in ticks:
        ax.text(3.55, t, f"{t * TR_S:.0f} s", fontsize=6, color=INK2, va="center")
    ax.axhline(79894, color=INK2, lw=1, ls="--")
    ax.text(-0.4, 79894 * 1.15, "no binning", fontsize=6.5, color=INK2)
    ax.set_xlim(-0.5, 3.5)
    ax.set_xticks(range(4)); ax.set_xticklabels(["Clean", "50%", "75%", "95%"], fontsize=7)
    ax.set_ylabel("readouts per gaze bin")
    ax.set_title("B", loc="left", fontsize=10, weight="bold", color=INK, x=-0.45)
    ax.text(1.5, 2 ** 17.6, "15 subjects x 4 bins\n(right: fixation time)", fontsize=6.5, color=INK2, ha="center")
    save(fig, "fig1_design_examples")


REGION_LABEL = {"brain": "Brain", "wm": "WM", "gm": "GM", "csf": "CSF", "volume": "Whole FoV",
                "eyes_mni": "MNI 'eye' block", "eyeball": "Eyeball", "vitreous": "Vitreous", "lens": "Lens",
                "nerve": "Optic nerve", "muscle": "Rectus muscles", "fat": "Orbital fat", "orbit": "Orbit"}


def subject_level():
    sl = pd.read_csv(AN / "subject_level.csv")
    return sl


def mean_sd(sl, source, method, region, metric, col="mean_bins"):
    x = sl[(sl.source == source) & (sl.method == method) & (sl.region == region) & (sl.metric == metric)]
    g = x.groupby("level")[col]
    return g.mean(), g.std(), g.count()


def fig2(regions=("brain", "wm", "eyeball", "lens", "nerve", "muscle")):
    """SSIM to the clean bin, per region, REGU vs GRID; band = no-binning image."""
    sl = subject_level()
    lv = ["50p", "75p", "95p"]
    fig, axs = plt.subplots(1, len(regions), figsize=(7.2, 2.2), sharey=True)
    for ax, reg in zip(axs, regions):
        for method, c, ls in (("REGU", C_REGU, "-"), ("GRID", C_GRID, "--")):
            m, sd, n = mean_sd(sl, "fixed", method, reg, "ssim")
            ax.errorbar(range(3), m[lv], yerr=sd[lv], color=c, ls=ls, marker="o", ms=4, lw=1.6,
                        capsize=2, elinewidth=1, label=method)
            wb, wsd = m["woBin"], sd["woBin"]
            if method == "REGU":
                ax.axhspan(wb - wsd, wb + wsd, color=INK2, alpha=0.12, lw=0)
                ax.axhline(wb, color=INK2, lw=1, ls=":")
        ax.set_title(REGION_LABEL[reg], color=INK)
        ax.set_xticks(range(3)); ax.set_xticklabels(["50%", "75%", "95%"])
        ax.set_xlim(-0.3, 2.3)
    axs[0].set_ylabel("SSIM to clean bin")
    axs[0].set_ylim(0.2, 1.0)
    axs[0].legend(loc="lower left", fontsize=7, handlelength=1.8)
    axs[-1].text(2.25, 0.21, "band: no binning\n(REGU, mean ± SD)", fontsize=6, color=INK2, ha="right", va="bottom")
    n = int(mean_sd(sl, "fixed", "REGU", regions[0], "ssim")[2].max())
    fig.suptitle(f"undersampling of the clean gaze bin  (n = {n} subjects, mean ± SD over subjects of the 4-bin mean)",
                 fontsize=7, color=INK2, y=1.02)
    fig.tight_layout(w_pad=0.4)
    save(fig, "fig2_ssim_regions")


def _knee_curve(d, k, b):
    """Population curve of the broken-stick fit: mean subject intercept + hinge."""
    h = np.minimum(np.log2(d.n_readouts.values) - k, 0)
    a = (d.value.values - b * h)
    return pd.Series(a).groupby(d.subject.values).mean().mean()


def fig3(panels=(("nerve_fat", "cnr", "Optic nerve–fat CNR"),
                 ("muscle_fat", "cnr", "Rectus–fat CNR"),
                 ("eyeball", "dcte_crop", "Eyeball DCTE (focus)"),
                 ("eyeball", "edge_width", "Eyeball edge width (mm)"),
                 ("vitreous", "cv", "Vitreous CV"))):
    """Readout budget: metric vs readouts per bin, all subjects x bins x levels.

    Broken-stick fit where it beats a straight line (F-test p < 0.05, knee inside
    the sampled range), otherwise the straight line, annotated per halving."""
    import runpy
    an = runpy.run_path(str(Path(__file__).with_name("07_analysis.py")))
    long = pd.read_csv(AN / "long.csv")
    kn = pd.read_csv(AN / "knee.csv")
    fig, axs = plt.subplots(1, len(panels), figsize=(8.4, 2.4))
    grid = np.arange(8.0, 14.01, 0.05)
    xs = np.linspace(7, 14.2, 200)
    for ax, (reg, met, title) in zip(axs, panels):
        for method, c, ls in (("REGU", INK, "-"), ("GRID", C_GRID, "--")):
            d = long[(long.source == "fixed") & (long.method == method) & (long.region == reg)
                     & (long.metric == met) & (long.level != "woBin")].dropna(subset=["value"])
            if method == "REGU":
                for lvl in ["full", "50p", "75p", "95p"]:
                    q = d[d.level == lvl]
                    ax.scatter(q.n_readouts, q.value, s=5, color=LEVEL_RAMP[lvl], alpha=0.75, lw=0, zorder=2)
            row = kn[(kn.source == "fixed") & (kn.method == method) & (kn.region == reg) & (kn.metric == met)].iloc[0]
            l2n, y, subj = np.log2(d.n_readouts.values), d.value.values, d.subject.values
            if row.p_vs_linear < 0.05 and not row.at_grid_edge:
                k, b, *_ = an["fit_knee"](l2n, y, subj, grid)
                a = _knee_curve(d, k, b)
                ax.plot(2 ** xs, a + b * np.minimum(xs - k, 0), color=c, ls=ls, lw=1.6, zorder=3)
                if method == "REGU":
                    ax.axvspan(row.knee_ci_lo, row.knee_ci_hi, color=C_REGU, alpha=0.12, lw=0, zorder=1)
                    ax.axvline(row.knee_readouts, color=C_REGU, lw=1, zorder=1)
                    ax.text(0.97, 0.05, f"REGU threshold\n{row.knee_readouts:,.0f} ({row.knee_seconds:.0f} s)",
                            transform=ax.transAxes, fontsize=6, color=C_REGU, ha="right")
            else:
                codes = pd.factorize(subj)[0]
                cen = lambda v: v - (np.bincount(codes, weights=v) / np.bincount(codes))[codes]
                b = (cen(l2n) @ cen(y)) / (cen(l2n) @ cen(l2n))
                a = np.mean(y) - b * np.mean(l2n)
                ax.plot(2 ** xs, a + b * xs, color=c, ls=ls, lw=1.6, zorder=3)
                if method == "REGU":
                    ax.text(0.97, 0.05 if b < 0 else 0.88, f"{-b:+.2g} per halving", transform=ax.transAxes,
                            fontsize=6, color=INK, ha="right")
        ax.set_xscale("log", base=2)
        ax.set_xticks([256, 1024, 4096, 16384])
        ax.set_xticklabels(["256\n2 s", "1k\n8 s", "4k\n33 s", "16k\n131 s"], fontsize=6.5)
        ax.minorticks_off()
        ax.set_title(title, fontsize=7.5, color=INK)
    handles = [Line2D([], [], ls="", marker="o", ms=3.5, color=LEVEL_RAMP[l], label=f"REGU {LABEL[l]}")
               for l in ["full", "50p", "75p", "95p"]]
    handles += [Line2D([], [], color=INK, lw=1.6, label="REGU fit"),
                Line2D([], [], color=C_GRID, lw=1.6, ls="--", label="GRID fit")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=6.5, bbox_to_anchor=(0.5, -0.06),
               handletextpad=0.3, columnspacing=1.0)
    fig.text(0.01, 0.0, "readouts per gaze bin\n(fixation time)", fontsize=6.5, color=INK2, va="bottom")
    fig.tight_layout(w_pad=0.6, rect=(0, 0.06, 1, 1))
    save(fig, "fig3_readout_budget")


MRIQC_SHOW = [("cjv", "CJV", -1), ("cnr", "CNR", 1), ("snr_wm", "SNR (WM)", 1), ("fber", "FBER", 1),
              ("efc", "EFC", -1), ("fwhm_avg", "FWHM (mm)", -1), ("inu_med", "INU", -1), ("wm2max", "WM2MAX", 1)]


def fig4():
    """MRIQC whole-image IQMs (row 1) and MReyeQC eye IQMs (row 2) vs level, REGU vs GRID."""
    sl = subject_level()
    panels = [("mriqc", "mriqc", "fwhm_avg", "FWHM (mm) ↓"), ("mriqc", "mriqc", "snrd_total", "SNR Dietrich ↑"),
              ("mriqc", "mriqc", "efc", "EFC ↓"), ("mriqc", "mriqc", "cjv", "CJV ↓"),
              ("mriqc", "mriqc", "cnr", "CNR ↑"),
              ("mreyeqc", "mreyeqc", "seg_snr_GLOBE", "Globe SNR ↑"), ("mreyeqc", "mreyeqc", "seg_snr_LENS", "Lens SNR ↑"),
              ("mreyeqc", "mreyeqc", "seg_snr_MUSCLE", "Rectus SNR ↑"), ("mreyeqc", "mreyeqc", "seg_cnr", "Lens–globe CNR ↑"),
              ("mreyeqc_model", "mreyeqc", "rating", "MReyeQC rating ↑")]
    fig, axs = plt.subplots(2, 5, figsize=(7.2, 3.3))
    for ax, (src, reg, met, lab) in zip(axs.flat, panels):
        for method, c, ls in (("REGU", C_REGU, "-"), ("GRID", C_GRID, "--")):
            m, sd, n = mean_sd(sl, src, method, reg, met)
            if not len(m):
                continue
            lv = [l for l in LEVELS if l in m.index]
            ax.errorbar([LEVELS.index(l) for l in lv], m[lv], yerr=sd[lv], color=c, ls=ls, marker="o", ms=3,
                        lw=1.4, capsize=1.5, elinewidth=0.8, label=method)
        ax.set_title(lab, fontsize=7, color=INK)
        ax.set_xticks(range(5)); ax.set_xticklabels(["NB", "C", "50", "75", "95"], fontsize=6)
        ax.tick_params(axis="y", labelsize=6)
    axs[0, 0].legend(fontsize=6, loc="best")
    axs[0, 0].set_ylabel("MRIQC (brain)", fontsize=7)
    axs[1, 0].set_ylabel("MReyeQC (eye)", fontsize=7)
    fig.text(0.5, -0.01, "NB no binning · C clean bin · 50/75/95 % of the clean readouts removed · mean ± SD over 15 subjects",
             ha="center", fontsize=6.5, color=INK2)
    fig.tight_layout(h_pad=0.8, w_pad=0.4)
    save(fig, "fig4_iqms")


def fig5():
    """Summary table: REGU, % change vs clean (mean over subjects), equivalence and significance."""
    lt = pd.read_csv(AN / "level_tests.csv")
    rows = [("fixed", "wm_gm", "cnr", "Brain  WM–GM CNR"), ("fixed", "wm", "cv", "Brain  WM CV"), ("fixed", "brain", "dcte_crop", "Brain  DCTE"),
            ("mriqc", "mriqc", "fwhm_avg", "MRIQC  FWHM"), ("mriqc", "mriqc", "snrd_total", "MRIQC  SNR Dietrich"),
            ("mriqc", "mriqc", "efc", "MRIQC  EFC"), ("mriqc", "mriqc", "cjv", "MRIQC  CJV"),
            ("mriqc", "mriqc", "cnr", "MRIQC  CNR"), ("mriqc", "mriqc", "snr_wm", "MRIQC  SNR WM"),
            ("fixed", "nerve_fat", "cnr", "Eye  optic nerve–fat CNR"), ("fixed", "muscle_fat", "cnr", "Eye  rectus–fat CNR"),
            ("fixed", "lens_vitreous", "cnr", "Eye  lens–vitreous CNR"), ("fixed", "vitreous", "cv", "Eye  vitreous CV"),
            ("fixed", "eyeball", "edge_width", "Eye  eyeball edge width"),
            ("fixed", "eyeball", "dcte_crop", "Eye  eyeball DCTE"),
            ("mreyeqc", "mreyeqc", "seg_snr_GLOBE", "MReyeQC  globe SNR"),
            ("mreyeqc", "mreyeqc", "seg_snr_LENS", "MReyeQC  lens SNR"),
            ("mreyeqc", "mreyeqc", "seg_snr_MUSCLE", "MReyeQC  rectus SNR")]
    cells, labels = [], []
    for src, reg, met, lab in rows:
        r = lt[(lt.source == src) & (lt.method == "REGU") & (lt.region == reg) & (lt.metric == met)].set_index("level")
        if r.empty:
            continue
        line = [(f"{r.loc['full', 'mean']:,.0f}" if abs(r.loc['full', 'mean']) >= 100 else f"{r.loc['full', 'mean']:.3g}") if "full" in r.index else "–"]
        for lvl in ["50p", "75p", "95p"]:
            if lvl not in r.index:
                line.append("–"); continue
            q = r.loc[lvl]
            tag = ("≈" if q["p_tost5"] < 0.05 else "") + ("*" if q["p_holm"] < 0.05 else "")
            line.append(f"{q['pct_change']:+.1f}% {tag}")
        cells.append(line); labels.append(lab)
    fig, ax = plt.subplots(figsize=(5.2, 0.24 * len(cells) + 0.6))
    ax.axis("off")
    # explicit bbox keeps the table inside the axes, below the title (loc="center" + scale overflowed it)
    t = ax.table(cellText=cells, rowLabels=labels, colLabels=["Clean", "50 %", "75 %", "95 %"],
                 cellLoc="center", rowLoc="left", bbox=[0.3, 0.0, 0.7, 1.0])
    t.auto_set_font_size(False); t.set_fontsize(7)
    for (i, j), c in t.get_celld().items():
        c.set_edgecolor(GRIDC)
        if i == 0:
            c.set_text_props(weight="bold", color=INK)
    ax.set_title("REGU: % change of the subject mean vs the clean bin (the reference)", fontsize=7, color=INK2, pad=4)
    fig.text(0.5, 0.02, "≈ change confirmed smaller than ±5 % (equivalence test, p<0.05)    "
             "* change confirmed not zero (Wilcoxon, Holm p<0.05)\n"
             "≈* small but consistent across subjects    no mark: neither could be shown",
             ha="center", va="top", fontsize=6.5, color=INK2, linespacing=1.5)
    save(fig, "fig5_summary_table")


def preview(sub="sub-013", pos="up"):
    """Preview figure: no caption, no small text, legible as a thumbnail.

    One clean gaze bin and its 50/75/95 % undersampled versions (REGU), axial
    slice through the optic nerves so globe, lens, recti and nerve are all in."""
    rc = pd.read_csv(ROOT / "readout_counts.csv")
    lv = {"clean": "full", "clean_0.50": "50p", "clean_0.75": "75p", "clean_0.95": "95p"}
    rc["level"] = rc.lvl.map(lv)
    pidx = {"up": 0, "down": 1, "left": 2, "right": 3}[pos]
    n_here = rc[(rc["sub"] == int(sub[-3:])) & (rc.region == pidx)].set_index("level").used
    L = np.asanyarray(nib.load(ROOT / "masks" / sub / "aeye_labels.nii.gz").dataobj)
    z = int(round(ndi.center_of_mass(L == 3)[2]))
    xs = np.nonzero(L.any(axis=(1, 2)))[0]
    ys = np.nonzero(np.isin(L, [1, 2, 3]).any(axis=(0, 2)))[0]
    sx, sy = slice(xs.min() - 4, xs.max() + 5), slice(ys.min() - 3, ys.max() + 6)
    titles = {"full": "Clean", "50p": "−50 %", "75p": "−75 %", "95p": "−95 %"}
    fig, axs = plt.subplots(1, 4, figsize=(9, 3.5))
    for ax, lvl in zip(axs, ["full", "50p", "75p", "95p"]):
        v = vol(sub, "REGU", pos, lvl)[sx, sy, z]
        ax.imshow(v.T, cmap="gray", origin="lower", vmin=0, vmax=np.percentile(v, 99.5), interpolation="nearest")
        ax.axis("off")
        ax.set_title(titles[lvl], fontsize=17, color=INK, pad=6, weight="bold")
        n = int(n_here[lvl])
        ax.text(0.5, -0.05, f"{n:,} readouts", transform=ax.transAxes, ha="center",
                va="top", fontsize=14, color=INK, weight="bold")
        ax.text(0.5, -0.27, f"{n * TR_S:.0f} s fixation", transform=ax.transAxes, ha="center",
                va="top", fontsize=14, color=C_REGU, weight="bold")
    fig.tight_layout(w_pad=0.4)
    save(fig, "preview")


FIGURES = {"preview": preview, "fig1": fig1, "fig2": fig2, "fig3": fig3, "fig4": fig4, "fig5": fig5}

if __name__ == "__main__":
    for name in (sys.argv[1:] or FIGURES):
        FIGURES[name]()
