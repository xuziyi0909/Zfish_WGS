"""Genome-wide figures for the WGS report (PNG for slides, PDF for papers).

manhattan        per-variant -log10 p that the mutant carries more alt reads than the control
variant_density  high-confidence mutant-specific variants per Mb
homozygosity     fraction of homozygous variant calls per Mb, per sample
homozygosity_runs share of each 1 Mb window in runs of homozygosity (bcftools roh), per sample
coverage         normalised read depth and mutant/control depth ratio (large deletions/duplications)
candidate_map    chromosomes with priority 1-2 variants and structural variants marked
summary          variant categories, effect types, allele fractions, SV sizes
"""

import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import hypergeom

sys.stderr = open(snakemake.log[0], "w")
inp, out = snakemake.input, snakemake.output
controls = list(snakemake.params.controls)
mutants = list(snakemake.params.mutants)
# categories counted as high-confidence candidates (mutant_enriched comes from sibling mode)
HIGH_CONF = ["mutant_specific", "mutant_enriched"]
outdir = snakemake.params.outdir

# ---------------------------------------------------------------- style
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"
CHROM_COLORS = ("#1c5cab", "#86b6ef")      # alternating chromosomes (one hue, two steps)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
NEG, POS = "#2a78d6", "#e34948"            # diverging: control-biased / mutant-biased
plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
    "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "figure.dpi": 100, "savefig.dpi": 200, "savefig.bbox": "tight",
})


def save(fig, name):
    fig.savefig(f"{outdir}/{name}.png")
    fig.savefig(f"{outdir}/{name}.pdf")
    plt.close(fig)
    print("wrote", name, file=sys.stderr)


def note(ax, text):
    ax.text(0.5, 0.5, text, transform=ax.transAxes, ha="center", va="center", color=INK2)
    ax.set_axis_off()


# ---------------------------------------------------------------- genome layout
fai = pd.read_csv(inp.fai, sep="\t", header=None, usecols=[0, 1], names=["chrom", "length"],
                  dtype={"chrom": str})
main = fai[(fai["length"] >= 10_000_000) & (fai["chrom"] != "MT")].copy()
# natural order: 1, 2, ... 25 (the FASTA index may list 1, 10, 11, ...)
main["_key"] = main["chrom"].map(lambda c: (0, int(c), "") if c.isdigit() else (1, 0, c))
main = main.sort_values("_key").drop(columns="_key").reset_index(drop=True)
GAP = 8_000_000
main["offset"] = (main["length"] + GAP).cumsum().shift(fill_value=0)
OFFSET = dict(zip(main["chrom"], main["offset"]))
LENGTH = dict(zip(main["chrom"], main["length"]))
CHROMS = list(main["chrom"])
GENOME_END = main["offset"].iloc[-1] + main["length"].iloc[-1]


def genome_x(chrom, pos):
    return chrom.map(OFFSET) + pos


def chrom_axis(ax):
    mids = [OFFSET[c] + LENGTH[c] / 2 for c in CHROMS]
    ax.set_xticks(mids)
    ax.set_xticklabels(CHROMS, fontsize=8)
    ax.set_xlim(-GAP, GENOME_END + GAP)
    ax.set_xlabel("chromosome")
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def chrom_color(chrom):
    idx = {c: i % 2 for i, c in enumerate(CHROMS)}
    return np.array(CHROM_COLORS)[chrom.map(idx).to_numpy()]


# ---------------------------------------------------------------- load variants
samples = mutants + controls
usecols = ["category", "chrom", "pos", "type", "impact"] + [
    f"{s}_{k}" for s in samples for k in ("genotype", "ref_reads", "alt_reads")]
dtypes = {"chrom": str, "category": "category", "type": "category", "impact": "category",
          **{f"{s}_genotype": "category" for s in samples},
          **{f"{s}_{k}": "int32" for s in samples for k in ("ref_reads", "alt_reads")}}
v = pd.read_csv(inp.variants, sep="\t", usecols=usecols, dtype=dtypes)
v = v[v["chrom"].isin(OFFSET)].reset_index(drop=True)
v["gx"] = genome_x(v["chrom"], v["pos"])
v["win"] = (v["pos"] // 1_000_000).astype("int32")
print(f"{len(v):,} variants on main chromosomes", file=sys.stderr)

mut = pd.read_csv(inp.mutant, sep="\t", dtype={"chrom": str},
                  usecols=["priority", "category", "chrom", "pos", "gene", "gene_id", "effect",
                           "impact", "hgvs_p"] + [f"{s}_alt_reads" for s in mutants]
                  + [f"{s}_ref_reads" for s in mutants])
mut = mut[mut["chrom"].isin(OFFSET)]
mut["gx"] = genome_x(mut["chrom"], mut["pos"])
top = mut[mut["priority"].isin([1, 2])].copy()
top["label"] = top["gene"].fillna(top["gene_id"]).astype(str)

m_alt = v[[f"{s}_alt_reads" for s in mutants]].sum(axis=1).to_numpy()
m_ref = v[[f"{s}_ref_reads" for s in mutants]].sum(axis=1).to_numpy()

def peak_region(df, col, frac=0.25):
    """Highest window, extended along its chromosome while values stay >= frac * peak."""
    if df[col].max() <= 0:
        return None, 0, 0, 0
    top_row = df.loc[df[col].idxmax()]
    chrom, peak = top_row["chrom"], top_row[col]
    d = df[df["chrom"] == chrom].set_index("win")[col]
    lo = hi = int(top_row["win"])
    while d.get(lo - 1, 0) >= frac * peak or d.get(lo - 2, 0) >= frac * peak:
        lo -= 1
    while d.get(hi + 1, 0) >= frac * peak or d.get(hi + 2, 0) >= frac * peak:
        hi += 1
    return chrom, lo, hi + 1, peak


# ================================================================ 1. Manhattan
fig, ax = plt.subplots(figsize=(13, 4.2))
if controls:
    c_alt = v[[f"{s}_alt_reads" for s in controls]].sum(axis=1).to_numpy()
    c_ref = v[[f"{s}_ref_reads" for s in controls]].sum(axis=1).to_numpy()
    n_mut, n_all = m_alt + m_ref, m_alt + m_ref + c_alt + c_ref
    k_alt = m_alt + c_alt
    keep = m_alt > 0
    # one-sided Fisher exact test: is the alt allele enriched in the mutant reads?
    p = np.ones(len(v))
    p[keep] = hypergeom.sf(m_alt[keep] - 1, n_all[keep], k_alt[keep], n_mut[keep])
    v["neglog10p"] = -np.log10(np.clip(p, 1e-300, 1))
    sig = v[v["neglog10p"] >= 1]
    rest = v[(v["neglog10p"] > 0) & (v["neglog10p"] < 1)].sample(frac=0.03, random_state=1) \
        if ((v["neglog10p"] > 0) & (v["neglog10p"] < 1)).any() else v.iloc[:0]
    pts = pd.concat([rest, sig])
    ax.scatter(pts["gx"], pts["neglog10p"], c=chrom_color(pts["chrom"]), s=2, linewidths=0,
               rasterized=True)
    if len(top):
        tp = top.merge(v[["chrom", "pos", "neglog10p"]], on=["chrom", "pos"], how="left")
        ax.scatter(tp["gx"], tp["neglog10p"], s=28, color=ORANGE, edgecolor="white",
                   linewidth=0.8, zorder=3)
        for _, r in tp[tp["priority"] == 1].iterrows():
            ax.annotate(r["label"], (r["gx"], r["neglog10p"]), xytext=(6, 14),
                        textcoords="offset points", ha="left", fontsize=8, color=INK,
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=GRID, lw=0.6),
                        arrowprops=dict(arrowstyle="-", color=INK2, lw=0.6))
    # shade the region where homozygous mutant-specific variants concentrate
    hc_hom = v[v["category"].isin(HIGH_CONF)
               & v[[f"{s}_genotype" for s in mutants]].eq("homozygous").any(axis=1)]
    _w = hc_hom.groupby(["chrom", "win"], observed=True).size().rename("n").reset_index()
    _c, _s, _e, _pk = peak_region(_w, "n")
    if _pk:
        ax.axvspan(OFFSET[_c] + _s * 1e6, OFFSET[_c] + _e * 1e6, color="#fbd5c4", zorder=0, lw=0,
                   label=f"homozygosity peak (chr{_c} {_s}–{_e} Mb)")
    chrom_axis(ax)
    ax.set_ylabel("−log10 p  (mutant alt-allele enrichment)")
    ax.set_ylim(bottom=0)
    ax.set_title("Manhattan plot: variants enriched in the mutant")
    handles = [Line2D([], [], marker="o", ls="", color=CHROM_COLORS[0], ms=4, label="variant"),
               Line2D([], [], marker="o", ls="", color=ORANGE, mec="white", ms=7,
                      label="priority 1–2 candidate")]
    handles += [h for h in ax.get_legend_handles_labels()[0] if h.get_label().startswith("homozygosity")]
    ax.legend(handles=handles,
              loc="upper right", frameon=False, fontsize=8)
    ax.text(0, -0.2, "One-sided Fisher exact test of alt vs ref read counts, mutant vs control. "
            "The highest reachable score depends on read depth, so low coverage limits the peak height.",
            transform=ax.transAxes, fontsize=8, color=INK2)
else:
    note(ax, "Manhattan plot needs at least one control sample.")
save(fig, "manhattan")

# ================================================================ 2. variant density (1 Mb)
win = v.groupby(["chrom", "win"], observed=True).size().rename("all_variants").to_frame()
hc = v["category"].isin(HIGH_CONF)
hom_m = v[[f"{s}_genotype" for s in mutants]].eq("homozygous").any(axis=1)
win["mutant_specific"] = v[hc].groupby(["chrom", "win"], observed=True).size()
win["mutant_specific_homozygous"] = v[hc & hom_m].groupby(["chrom", "win"], observed=True).size()
win = win.fillna(0).astype(int).reset_index()
PEAK = peak_region(win, "mutant_specific_homozygous")

fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True)
for ax, col, title in [
    (axes[0], "mutant_specific", "High-confidence mutant-specific variants per Mb"),
    (axes[1], "mutant_specific_homozygous", "…of which homozygous in the mutant"),
]:
    x = genome_x(win["chrom"], win["win"] * 1_000_000 + 500_000)
    ax.bar(x, win[col], width=900_000, color=chrom_color(win["chrom"]), linewidth=0)
    ax.set_ylabel("variants / Mb")
    ax.set_title(title)
    chrom_axis(ax)
    c0, s0, e0, peak = peak_region(win, col)
    if peak > 0:
        ax.annotate(f"peak: chr{c0} {s0}–{e0} Mb", (OFFSET[c0] + s0 * 1e6, peak),
                    xytext=(-8, -4), textcoords="offset points", ha="right", va="top",
                    fontsize=8.5, color=INK, fontweight="bold")
axes[0].set_xlabel("")
fig.tight_layout()
save(fig, "variant_density")

# ================================================================ 3. homozygosity (1 Mb)
fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True,
                         gridspec_kw={"height_ratios": [1.3, 1]})
hz = {}
for s in samples:
    depth = v[f"{s}_ref_reads"] + v[f"{s}_alt_reads"]
    called = (depth >= 4) & (v[f"{s}_alt_reads"] > 0)
    is_hom = v[f"{s}_alt_reads"] == depth
    g = pd.DataFrame({"chrom": v["chrom"], "win": v["win"], "called": called, "hom": called & is_hom})
    agg = g.groupby(["chrom", "win"], observed=True)[["called", "hom"]].sum()
    hz[s] = (agg["hom"] / agg["called"].where(agg["called"] >= 20)).rename(s)
hzdf = pd.concat(hz.values(), axis=1).reset_index()
hzdf["x"] = genome_x(hzdf["chrom"], hzdf["win"] * 1_000_000 + 500_000)
for s in (controls + mutants)[:4]:          # mutants drawn last = on top
    color = ORANGE if s in mutants else BLUE
    for c in CHROMS:
        d = hzdf[hzdf["chrom"] == c]
        smooth = d[s].rolling(5, center=True, min_periods=2).mean()
        axes[0].plot(d["x"], smooth, color=color, lw=1.6, label=s if c == CHROMS[0] else None)
axes[0].set_ylim(0, 1.02)
axes[0].set_ylabel("fraction homozygous")
axes[0].set_title("Homozygosity along the genome (sites with ≥ 4 reads, 5 Mb running mean)")
axes[0].legend(frameon=False, fontsize=8, loc="upper right", ncol=len(samples[:3]))
chrom_axis(axes[0])
axes[0].set_xlabel("")
if controls:
    diff = hzdf[mutants].mean(axis=1) - hzdf[controls].mean(axis=1)
    axes[1].bar(hzdf["x"], diff, width=900_000, linewidth=0,
                color=np.where(diff >= 0, POS, NEG))
    axes[1].axhline(0, color=INK2, lw=0.8)
    axes[1].set_ylabel("mutant − control")
    axes[1].set_title("Difference: red = more homozygous in the mutant")
    chrom_axis(axes[1])
    hzdf["mutant_minus_control"] = diff
else:
    note(axes[1], "Difference track needs a control sample.")
fig.tight_layout()
save(fig, "homozygosity")

# ================================================================ 3b. runs of homozygosity (bcftools roh)
ROH_PEAK = None
roh = None
if inp.roh:
    runs = pd.read_csv(inp.roh, sep="\t", dtype={"chrom": str, "sample": str})
    runs = runs[runs["chrom"].isin(OFFSET) & runs["sample"].isin(samples)]
    roh = pd.DataFrame([(c, w) for c in CHROMS for w in range(-(-LENGTH[c] // 1_000_000))],
                       columns=["chrom", "win"])
    win_size = (np.minimum((roh["win"] + 1) * 1_000_000, roh["chrom"].map(LENGTH))
                - roh["win"] * 1_000_000)
    # base pairs of each run that fall into each 1 Mb window
    parts = [(r.sample, r.chrom, w, min(r.end, (w + 1) * 1_000_000) - max(r.start, w * 1_000_000))
             for r in runs.itertuples()
             for w in range(r.start // 1_000_000, (r.end - 1) // 1_000_000 + 1)]
    parts = pd.DataFrame(parts, columns=["sample", "chrom", "win", "bp"])
    for s in samples:
        bp = parts[parts["sample"] == s].groupby(["chrom", "win"])["bp"].sum().rename(s).reset_index()
        roh = roh.merge(bp, on=["chrom", "win"], how="left")
        roh[s] = (roh[s].fillna(0) / win_size).clip(0, 1)
    roh["x"] = genome_x(roh["chrom"], roh["win"] * 1_000_000 + 500_000)

    fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.3, 1]})
    for s in (controls + mutants)[:4]:          # mutants drawn last = on top
        color = ORANGE if s in mutants else BLUE
        for c in CHROMS:
            d = roh[roh["chrom"] == c]
            axes[0].plot(d["x"], d[s].rolling(3, center=True, min_periods=1).mean(), color=color, lw=1.4,
                         label=s if c == CHROMS[0] else None)
    axes[0].set_ylim(0, 1.02)
    axes[0].set_ylabel("share of window in runs")
    axes[0].set_title("Runs of homozygosity (bcftools roh): share of each 1 Mb window, 3 Mb running mean")
    axes[0].legend(frameon=False, fontsize=8, loc="upper right", ncol=len(samples[:3]))
    chrom_axis(axes[0])
    axes[0].set_xlabel("")
    if controls:
        roh["mutant_minus_control"] = roh[mutants].mean(axis=1) - roh[controls].mean(axis=1)
        diff = roh["mutant_minus_control"]
        axes[1].bar(roh["x"], diff, width=900_000, linewidth=0, color=np.where(diff >= 0, POS, NEG))
        axes[1].axhline(0, color=INK2, lw=0.8)
        axes[1].set_ylim(-1.05, 1.05)
        axes[1].set_ylabel("mutant − control")
        axes[1].set_title("Difference: red = the mutant is in runs where the control is not")
        chrom_axis(axes[1])
        # strongest block: the highest window, extended while neighbours stay above half of it
        c0, s0, e0, top_ = peak_region(roh, "mutant_minus_control", frac=0.5)
        if top_ >= 0.3:
            inside = roh[(roh["chrom"] == c0) & (roh["win"] >= s0) & (roh["win"] < e0)]
            ROH_PEAK = {"chrom": c0, "start_mb": int(s0), "end_mb": int(e0),
                        "max_difference": round(float(top_), 2),
                        "mutant_share": round(float(inside[mutants].mean(axis=1).mean()), 2),
                        "control_share": round(float(inside[controls].mean(axis=1).mean()), 2)}
            axes[1].axvspan(OFFSET[c0] + s0 * 1e6, OFFSET[c0] + e0 * 1e6, color="#fbd5c4", zorder=0, lw=0)
            axes[1].annotate(f"chr{c0} {s0}–{e0} Mb", (OFFSET[c0] + e0 * 1e6, top_), xytext=(6, -4),
                             textcoords="offset points", ha="left", va="top", fontsize=8.5,
                             color=INK, fontweight="bold")
    else:
        note(axes[1], "Difference track needs a control sample.")
    fig.tight_layout()
    save(fig, "homozygosity_runs")

win = win.merge(hzdf.drop(columns="x").rename(columns={s: f"{s}_fraction_homozygous" for s in samples}),
                on=["chrom", "win"], how="left")
if roh is not None:
    roh_cols = list(samples) + (["mutant_minus_control"] if controls else [])
    win = win.merge(roh[["chrom", "win"] + roh_cols].rename(
        columns={**{s: f"{s}_fraction_in_runs" for s in samples},
                 "mutant_minus_control": "runs_mutant_minus_control"}), on=["chrom", "win"], how="left")
win.insert(1, "start", win.pop("win") * 1_000_000)
win.insert(2, "end", win["start"] + 1_000_000)
win.round(3).to_csv(out.windows, sep="\t", index=False)

# ================================================================ 4. coverage / copy number
cov = pd.read_csv(inp.coverage, sep="\t", dtype={"chrom": str})
cov = cov[cov["chrom"].isin(OFFSET)]
size = (cov["end"] - cov["start"]).to_numpy()
cov["x"] = genome_x(cov["chrom"], (cov["start"] + cov["end"]) // 2)
norm = {}
for s in samples:
    depth = cov[s] / size
    norm[s] = depth / np.median(depth[depth > 0])
fig, axes = plt.subplots(2, 1, figsize=(13, 5.6), sharex=True)
for s in samples[:3]:
    color = ORANGE if s in mutants else BLUE
    axes[0].scatter(cov["x"], norm[s], s=3, color=color, linewidths=0, label=s, rasterized=True)
axes[0].set_ylim(0, 2.5)
axes[0].axhline(1, color=INK2, lw=0.8)
axes[0].set_ylabel("depth / genome median")
axes[0].set_title("Read depth in 500 kb windows")
axes[0].legend(frameon=False, fontsize=8, loc="upper right", markerscale=3, ncol=3)
chrom_axis(axes[0])
axes[0].set_xlabel("")
if controls:
    m = np.mean([norm[s] for s in mutants], axis=0)
    c = np.mean([norm[s] for s in controls], axis=0)
    ok = (m > 0) & (c > 0.2)
    ratio = pd.Series(np.where(ok, np.log2(np.where(ok, m, 1) / np.where(ok, c, 1)), np.nan))
    smooth = ratio.groupby(cov["chrom"].to_numpy()).transform(
        lambda r: r.rolling(5, center=True, min_periods=3).median())
    axes[1].scatter(cov["x"], ratio, s=3, c=chrom_color(cov["chrom"]), linewidths=0, rasterized=True)
    for ch in CHROMS:
        sel = (cov["chrom"] == ch).to_numpy()
        axes[1].plot(cov["x"].to_numpy()[sel], smooth.to_numpy()[sel], color=INK, lw=1)
    for y in (-1, 0, 0.585):
        axes[1].axhline(y, color=INK2 if y == 0 else GRID, lw=0.8, ls="-" if y == 0 else "--")
    axes[1].set_ylim(-2, 2)
    axes[1].set_ylabel("log2(mutant / control)")
    axes[1].set_title("Copy-number ratio  (−1 ≈ one copy lost, +0.58 ≈ one copy gained)")
    chrom_axis(axes[1])
    cov["log2_mutant_over_control"] = ratio.round(3).to_numpy()
else:
    note(axes[1], "Copy-number ratio needs a control sample.")
for s in samples:
    cov[f"{s}_mean_depth"] = (cov[s] / size).round(2)
cov.drop(columns=["x"] + samples).to_csv(out.cov_ratio, sep="\t", index=False)
fig.tight_layout()
save(fig, "coverage")

# ================================================================ 5. candidate map
sv = pd.DataFrame()
if inp.sv:
    try:
        sv = pd.read_csv(inp.sv, sep="\t", dtype={"chrom": str})
        sv = sv[sv["chrom"].isin(OFFSET)]
    except pd.errors.EmptyDataError:
        sv = pd.DataFrame()
fig, ax = plt.subplots(figsize=(11, 8))
for i, c in enumerate(CHROMS):
    ax.barh(i, LENGTH[c] / 1e6, height=0.5, color="#ecebe7", edgecolor=GRID, linewidth=0.6)
ypos = {c: i for i, c in enumerate(CHROMS)}
if PEAK[0] is not None:
    pc_, ps_, pe_, _ = PEAK
    ax.add_patch(plt.Rectangle((ps_, ypos[pc_] - 0.42), pe_ - ps_, 0.84, color="#fbd5c4", zorder=1,
                               label="peak of homozygous mutant-specific variants"))
if len(sv):
    gene_sv = sv[sv["genes_overlapped"].notna() & (sv["genes_overlapped"] != "")]
    ax.scatter(gene_sv["start"] / 1e6, gene_sv["chrom"].map(ypos) + 0.36, marker="v", s=22,
               color=INK2, label="structural variant overlapping a gene", zorder=3)
p2 = top[top["priority"] == 2]
p1 = top[top["priority"] == 1]
ax.vlines(p2["pos"] / 1e6, p2["chrom"].map(ypos) - 0.25, p2["chrom"].map(ypos) + 0.25,
          color=BLUE, lw=1.4, label="priority 2 variant", zorder=3)
ax.scatter(p1["pos"] / 1e6, p1["chrom"].map(ypos), s=70, color=ORANGE, edgecolor="white",
           linewidth=1, label="priority 1 variant", zorder=4)
for _, r in p1.iterrows():
    ax.annotate(r["label"], (r["pos"] / 1e6, ypos[r["chrom"]]), xytext=(6, -12),
                textcoords="offset points", fontsize=8.5, color=INK, fontweight="bold")
ax.set_yticks(range(len(CHROMS)))
ax.set_yticklabels([f"chr{c}" for c in CHROMS], fontsize=8.5)
ax.invert_yaxis()
ax.set_xlabel("position (Mb)")
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.set_title("Where are the candidates?")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
save(fig, "candidate_map")

# ================================================================ 6. summary panel
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
# a. categories
CAT_STYLE = [  # (category, label, colour); categories absent from the data are skipped
    ("mutant_specific", "mutant-specific (high conf.)", ORANGE),
    ("mutant_enriched", "mutant-enriched (sibling mode)", "#e8875a"),
    ("mutant_specific_low_confidence", "mutant-specific (low conf.)", "#f5b294"),
    ("mutant_heterozygous", "mutant, not homozygous", "#d9cfc4"),
    ("shared", "shared", BLUE),
    ("control_only", "control only", "#86b6ef"),
    ("low_evidence", "low evidence", "#c3c2b7"),
]
present = set(v["category"].unique())
CAT_STYLE = [c for c in CAT_STYLE if c[0] in present or c[0] in ("mutant_specific", "shared")]
cats, labels, colors = zip(*CAT_STYLE)
counts = v["category"].value_counts().reindex(list(cats), fill_value=0)
ax = axes[0, 0]
ax.barh(list(labels), counts.to_numpy() / 1e6, color=list(colors), height=0.6)
for i, n in enumerate(counts):
    ax.text(n / 1e6, i, f"  {n:,}", va="center", fontsize=8.5, color=INK)
ax.invert_yaxis()
ax.set_xlabel("variants (millions)")
ax.set_title("a. How variants were classified (chr 1–25)")
ax.set_xlim(0, counts.max() / 1e6 * 1.3)
# b. effect types among protein-changing mutant-specific
ax = axes[0, 1]
pc = mut[mut["impact"].isin(["HIGH", "MODERATE"]) & mut["category"].isin(HIGH_CONF)]
eff = pc["effect"].str.split(", ").str[0].value_counts().head(10)[::-1]
if len(eff):
    ax.barh([e.replace("_", " ") for e in eff.index], eff.to_numpy(), color=BLUE, height=0.6)
    for i, n in enumerate(eff):
        ax.text(n, i, f" {n:,}", va="center", fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ax.set_xlabel("variants (log scale)")
else:
    note(ax, "No protein-changing mutant-specific variants.")
ax.set_title("b. Protein-changing effects (high confidence)")
# c. alt allele fraction of mutant-specific variants in the mutant
ax = axes[1, 0]
hcm = mut[mut["category"].isin(HIGH_CONF)]
alt = hcm[[f"{s}_alt_reads" for s in mutants]].sum(axis=1)
dep = alt + hcm[[f"{s}_ref_reads" for s in mutants]].sum(axis=1)
frac = (alt / dep.where(dep >= 4)).dropna()
ax.hist(frac, bins=np.linspace(-0.025, 1.025, 22), color=ORANGE, edgecolor="white", linewidth=1)
ax.set_xlabel("fraction of mutant reads carrying the variant\n"
              "≈0.5 = one copy (heterozygous)   ≈1.0 = both copies (homozygous)")
ax.set_ylabel("variants")
ax.set_title("c. Allele fraction in the mutant (sites ≥ 4 reads)")
# d. SV sizes by type
ax = axes[1, 1]
if len(sv) and sv["size_bp"].notna().any():
    for t, color in zip(["DEL", "DUP", "INV", "INS"], [ORANGE, BLUE, AQUA, "#6250d6"]):
        s_ = sv.loc[(sv["type"] == t) & (sv["size_bp"] > 0), "size_bp"]
        if len(s_):
            ax.hist(np.log10(s_), bins=np.linspace(1, 8, 29), color=color, alpha=0.85,
                    label=f"{t} ({len(s_)})", edgecolor="white", linewidth=0.6)
    ax.set_xticks(range(1, 9))
    ax.set_xticklabels(["10 bp", "100", "1 kb", "10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb"],
                       fontsize=8)
    ax.set_xlabel("structural variant size")
    ax.set_ylabel("variants")
    ax.legend(frameon=False, fontsize=8.5)
else:
    note(ax, "No mutant-specific structural variants.")
ax.set_title("d. Structural variant sizes")
fig.tight_layout()
save(fig, "summary")

# ---------------------------------------------------------------- peak summary for the report

peak = {"found": False}
if PEAK[0] is not None:
    pc_, ps_, pe_, ph_ = PEAK
    sel = (win["chrom"] == pc_) & (win["start"] >= ps_ * 1e6) & (win["end"] <= pe_ * 1e6)
    in_peak = top[(top["chrom"] == pc_) & (top["pos"] >= ps_ * 1e6) & (top["pos"] < pe_ * 1e6)]
    all_in_peak = mut[(mut["priority"] > 0) & (mut["chrom"] == pc_)
                      & (mut["pos"] >= ps_ * 1e6) & (mut["pos"] < pe_ * 1e6)]
    all_in_peak.drop(columns=["gx"]).sort_values(["priority", "pos"]).to_csv(
        out.peak_candidates, sep="\t", index=False)
    genome_mean = win["mutant_specific_homozygous"].mean()
    peak = {
        "found": True, "chrom": pc_, "start_mb": int(ps_), "end_mb": int(pe_),
        "peak_per_mb": int(ph_), "genome_mean_per_mb": round(float(genome_mean), 1),
        "homozygous_variants_in_peak": int(win.loc[sel, "mutant_specific_homozygous"].sum()),
        "homozygous_variants_total": int(win["mutant_specific_homozygous"].sum()),
        "candidates": in_peak.sort_values(["priority", "pos"])[
            ["priority", "chrom", "pos", "label", "gene_id", "effect", "hgvs_p"]]
            .fillna("").to_dict("records"),
    }
if not peak["found"]:
    pd.DataFrame(columns=["priority", "chrom", "pos", "gene"]).to_csv(out.peak_candidates, sep="\t", index=False)
peak["runs_of_homozygosity"] = ROH_PEAK
with open(out.peak, "w") as fh:
    json.dump(peak, fh, indent=1, default=str)
