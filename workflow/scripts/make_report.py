"""Build a single self-contained HTML report written for non-specialists."""

import base64
import datetime
import html
import json
import re
import sys

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
inp = snakemake.input
samples = snakemake.params.samples
cfg = snakemake.params.config
names = [s["sample"] for s in samples]
groups = {s["sample"]: s["group"] for s in samples}
paired = {s["sample"]: bool(s.get("fastq_2")) for s in samples}
esc = html.escape


# ------------------------------------------------------------------ QC metrics
def parse_flagstat(path):
    for line in open(path):
        if "primary mapped" in line:
            m = re.search(r"\(([\d.]+)%", line)
            return float(m.group(1)) if m else None
    return None


def parse_markdup(path):
    vals = {}
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            vals[k.strip()] = v.strip()
    examined = float(vals.get("EXAMINED", 0) or 0)
    dups = float(vals.get("DUPLICATE TOTAL", 0) or 0)
    return 100 * dups / examined if examined else None


def parse_coverage(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]
    df["length"] = df["endpos"] - df["startpos"] + 1
    total = df["length"].sum()
    return {
        "mean_depth": (df["meandepth"] * df["length"]).sum() / total,
        "breadth": 100 * df["covbases"].sum() / total,
        "per_chrom": df,
    }


qc = []
cov_tables = {}
for i, s in enumerate(names):
    fp = json.load(open(inp.fastp[i]))
    cov = parse_coverage(inp.coverage[i])
    cov_tables[s] = cov["per_chrom"]
    qc.append({
        "sample": s,
        "group": groups[s],
        "reads": "paired-end" if paired[s] else "single-end",
        "raw_reads": fp["summary"]["before_filtering"]["total_reads"],
        "kept_pct": 100 * fp["summary"]["after_filtering"]["total_reads"]
                    / max(fp["summary"]["before_filtering"]["total_reads"], 1),
        "q30": 100 * fp["summary"]["after_filtering"]["q30_rate"],
        "mapped_pct": parse_flagstat(inp.flagstat[i]),
        "dup_pct": parse_markdup(inp.markdup[i]),
        "mean_depth": cov["mean_depth"],
        "breadth": cov["breadth"],
    })


def light(value, good, ok, higher_is_better=True):
    if value is None:
        return "na", "–"
    if not higher_is_better:
        value_cmp, good, ok = -value, -good, -ok
    else:
        value_cmp = value
    if value_cmp >= good:
        return "good", "Good"
    if value_cmp >= ok:
        return "warn", "Caution"
    return "bad", "Poor"


def status(cls, text, value):
    icon = {"good": "✓", "warn": "!", "bad": "✕", "na": "–"}[cls]
    return (f'<td class="num"><span class="val">{value}</span>'
            f'<span class="pill {cls}"><span aria-hidden="true">{icon}</span> {text}</span></td>')


qc_rows = []
for r in qc:
    qc_rows.append(
        "<tr>"
        f"<td><b>{esc(r['sample'])}</b></td><td>{esc(r['group'])}</td><td>{r['reads']}</td>"
        f'<td class="num">{r["raw_reads"]/1e6:,.1f} M</td>'
        + status(*light(r["kept_pct"], 90, 75), f"{r['kept_pct']:.1f}%")
        + status(*light(r["q30"], 85, 75), f"{r['q30']:.1f}%")
        + status(*light(r["mapped_pct"], 95, 85), f"{r['mapped_pct']:.1f}%")
        + status(*light(r["dup_pct"], 15, 30, higher_is_better=False), f"{r['dup_pct']:.1f}%")
        + status(*light(r["mean_depth"], 15, 8), f"{r['mean_depth']:.1f}×")
        + status(*light(r["breadth"], 95, 85), f"{r['breadth']:.1f}%")
        + "</tr>"
    )

# ------------------------------------------------------------------ coverage chart (inline SVG)
SERIES = ["var(--series-1)", "var(--series-2)", "var(--series-3)"]
main_chroms = [str(c) for c in range(1, 26)]


def coverage_svg():
    W, H, L, R, T, B = 760, 260, 44, 90, 14, 34
    vals = {s: cov_tables[s].set_index("rname")["meandepth"] for s in names}
    chroms = [c for c in main_chroms if all(c in vals[s].index for s in names)]
    if not chroms:
        return "<p>No main chromosomes found.</p>"
    ymax = max(max(vals[s][c] for c in chroms) for s in names) * 1.1 or 1
    step = (W - L - R) / max(len(chroms) - 1, 1)
    x = lambda i: L + i * step
    y = lambda v: T + (H - T - B) * (1 - v / ymax)
    out = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-label="Mean coverage per chromosome">']
    nice = [t for t in (1, 2, 5, 10, 20, 50, 100) if ymax / t <= 6][0]
    t = 0
    while t <= ymax:
        out.append(f'<line class="grid" x1="{L}" x2="{W-R}" y1="{y(t):.1f}" y2="{y(t):.1f}"/>'
                   f'<text class="tick" x="{L-6}" y="{y(t)+4:.1f}" text-anchor="end">{t}×</text>')
        t += nice
    for i, c in enumerate(chroms):
        out.append(f'<text class="tick" x="{x(i):.1f}" y="{H-B+16}" text-anchor="middle">{c}</text>')
    out.append(f'<text class="tick" x="{(L+W-R)/2}" y="{H-2}" text-anchor="middle">chromosome</text>')
    show = names[:3]
    for k, s in enumerate(show):
        pts = " ".join(f"{x(i):.1f},{y(vals[s][c]):.1f}" for i, c in enumerate(chroms))
        out.append(f'<polyline fill="none" stroke="{SERIES[k]}" stroke-width="2" points="{pts}"/>')
        for i, c in enumerate(chroms):
            out.append(f'<circle cx="{x(i):.1f}" cy="{y(vals[s][c]):.1f}" r="4" fill="{SERIES[k]}" '
                       f'class="dot"><title>{esc(s)} · chr{c}: {vals[s][c]:.1f}×</title></circle>')
        last = vals[s][chroms[-1]]
        out.append(f'<text class="lbl" x="{x(len(chroms)-1)+8:.1f}" y="{y(last)+4 + 12*k:.1f}">{esc(s)}</text>')
    out.append("</svg>")
    legend = "".join(f'<span class="key"><i style="background:{SERIES[k]}"></i>{esc(s)}</span>'
                     for k, s in enumerate(show))
    return f'<div class="legend">{legend}</div>' + "".join(out)


# ------------------------------------------------------------------ variants
summary = json.load(open(inp.summary))
cats = summary["categories"]
mode = summary.get("mode", "absent")
HIGH_CONF = ("mutant_specific", "mutant_enriched")
sibling = mode == "sibling"
cat_text = [
    ("mutant_specific", "Mutant-specific (high confidence)",
     ("In all mutant reads; " if sibling else "In the mutant; ")
     + "control is well covered and shows no trace of it."),
    ("mutant_enriched", "Mutant-enriched (high confidence)",
     "In all mutant reads and in some control reads – what a recessive mutation looks like when the "
     "control siblings carry one copy."),
    ("mutant_specific_low_confidence", "Mutant-specific (low confidence)",
     "In the mutant; control had too few reads to be sure."),
    ("mutant_heterozygous", "Mutant, not homozygous",
     "In some but not all mutant reads – cannot be a recessive mutation (sibling mode)."),
    ("shared", "Shared", "Present in both mutant and control – background variation, not the mutation."),
    ("control_only", "Control only", "Seen only in the control."),
    ("low_evidence", "Low evidence", "Too few supporting reads to classify."),
]
cat_text = [c for c in cat_text if c[0] in cats or c[0] in ("mutant_specific", "shared")]
cat_max = max([cats.get(k, 0) for k, *_ in cat_text] + [1])
cat_rows = "".join(
    f'<tr><td><b>{t}</b><div class="sub">{d}</div></td><td class="num">{cats.get(k,0):,}</td>'
    f'<td class="barcell"><div class="bar" style="width:{100*cats.get(k,0)/cat_max:.1f}%" '
    f'title="{t}: {cats.get(k,0):,}"></div></td>'
    f'<td class="num">{summary["impact_by_category"].get(k,{}).get("HIGH",0):,}</td>'
    f'<td class="num">{summary["impact_by_category"].get(k,{}).get("MODERATE",0):,}</td></tr>'
    for k, t, d in cat_text)

genes = pd.read_csv(inp.genes, sep="\t", dtype=str).fillna("")
_pk = json.load(open(inp.peak))


def in_peak(chrom, pos):
    return (_pk.get("found") and str(chrom) == str(_pk["chrom"])
            and _pk["start_mb"] * 1e6 <= int(pos) < _pk["end_mb"] * 1e6)


def peak_tag(chrom, pos):
    return ' <span class="pill p2">in peak</span>' if in_peak(chrom, pos) else ""
mutants = summary["mutants"]


def prio_badge(p):
    p = str(p)
    cls = {"1": "p1", "2": "p2"}.get(p, "na")
    return f'<span class="pill {cls}">{esc(p)}</span>'


def gene_link(name, gid):
    label = esc(name or gid)
    if gid.startswith("ENSDARG"):
        return (f'<a href="https://www.ensembl.org/Danio_rerio/Gene/Summary?g={gid}" target="_blank">{label}</a>'
                f' · <a class="sub" href="https://zfin.org/search?q={esc(name or gid)}" target="_blank">ZFIN</a>')
    return label


gene_rows = "".join(
    f"<tr><td class='num'>{prio_badge(g.best_priority)}</td><td>{gene_link(g.gene, g.gene_id)}{peak_tag(g.chrom, g.first_pos)}</td><td>{g.chrom}</td>"
    f'<td class="num">{g.high_impact_variants}</td><td class="num">{g.moderate_impact_variants}</td>'
    f'<td class="num">{g.high_confidence_variants}</td><td class="num">{g.homozygous_in_mutant}</td>'
    f"<td>{esc(g.effects)}</td></tr>"
    for g in genes.head(50).itertuples())

var = pd.read_csv(inp.variants, sep="\t", dtype=str, nrows=200000).fillna("")
top = var[var["priority"].isin(["1", "2"])].head(100)
CHECK_SHORT = {"matches_zebrafish_refseq": "matches zebrafish RefSeq",
               "matches_vertebrate_consensus": "usual vertebrate residue",
               "alt_seen_in_orthologs": "seen in other vertebrates",
               "conserved": "conserved in vertebrates", "variable": "variable in vertebrates"}


def checks_cell(v):
    bits = []
    if summary.get("known_variants_checked"):
        bits.append(f"known {esc(v.known_id)}" if v.known_id else "novel")
    if CHECK_SHORT.get(v.allele_check):
        bits.append(esc(CHECK_SHORT[v.allele_check]))
    return "<br>".join(bits)


def loc_link(chrom, pos):
    p = int(pos)
    return (f'<a href="https://www.ensembl.org/Danio_rerio/Location/View?r={chrom}:{max(p-50,1)}-{p+50}" '
            f'target="_blank">{chrom}:{p:,}</a>')


var_rows = "".join(
    f"<tr><td class='num'>{prio_badge(v.priority)}</td><td>{loc_link(v.chrom, v.pos)}{peak_tag(v.chrom, v.pos)}</td><td>{gene_link(v.gene, v.gene_id)}</td>"
    f"<td>{esc(v.effect)}</td><td><code>{esc(v.hgvs_p or v.hgvs_c)}</code></td><td class='sub'>{checks_cell(v)}</td>"
    + "".join(f"<td>{esc(getattr(v, s + '_genotype'))} "
              f"({getattr(v, s + '_alt_reads')}/{int(getattr(v, s + '_alt_reads') or 0) + int(getattr(v, s + '_ref_reads') or 0)})</td>"
              for s in names)
    + "</tr>"
    for v in top.itertuples())
var_head = "".join(f"<th>{esc(s)}<div class='sub'>genotype (alt/total reads)</div></th>" for s in names)

# ------------------------------------------------------------------ sanity checks on candidates
checks_html = ""
pch = summary.get("protein_changing_high_confidence", {})
if pch.get("total") and (summary.get("known_variants_checked") or summary.get("protein_check")):
    ac = pch.get("allele_check", {})
    items = []
    if summary.get("known_variants_checked"):
        items.append(f"<li><b>{pch['known']:,}</b> of {pch['total']:,} are <b>known</b> zebrafish variants "
                     "(Ensembl/dbSNP). These are strain differences; a new mutation is almost always novel.</li>")
    if summary.get("protein_check"):
        items.append(f"<li><b>{ac.get('matches_zebrafish_refseq', 0):,}</b> change the amino acid to the one found in a "
                     "published zebrafish protein (RefSeq, from cDNA): the reference genome probably carries a rare "
                     "allele or an error there, and the mutant is normal.</li>")
        items.append(f"<li><b>{ac.get('matches_vertebrate_consensus', 0):,}</b> change it to the residue most other "
                     "vertebrates have: probably harmless.</li>")
        items.append(f"<li><b>{ac.get('conserved', 0):,}</b> hit an amino acid that is identical in all compared "
                     "vertebrates: these are the most likely to be damaging.</li>")
    demote = cfg["comparison"].get("demote_flagged", True)
    checks_html = ('<div class="card" style="margin-top:12px"><b>Checks on the high-confidence protein-changing '
                   f'variants ({pch["total"]:,})</b><ul>{"".join(items)}</ul><div class="sub">'
                   + ("Known variants and the first two groups are moved to priority 5." if demote else
                      "Flagged variants keep their priority (demote_flagged: false); see the 'warnings' column.")
                   + "</div></div>")

# ------------------------------------------------------------------ structural variants
sv_html = "<p class='sub'>Structural-variant calling was switched off in config/config.yaml.</p>"
n_sv = None
if inp.sv:
    sv = pd.read_csv(inp.sv, sep="\t", dtype=str).fillna("") if open(inp.sv).read().strip() else pd.DataFrame()
    n_sv = len(sv)
    if n_sv:
        sv["_g"] = sv["genes_overlapped"] != ""
        sv = sv.sort_values("_g", ascending=False)
        sv_html = ("<table><thead><tr><th>Location</th><th>Type</th><th class='num'>Size (bp)</th>"
                   "<th class='num'>Split-read support</th><th>Genes overlapped</th></tr></thead><tbody>"
                   + "".join(f"<tr><td>{loc_link(r.chrom, r.start)}</td><td>{esc(r.type)}</td>"
                             f"<td class='num'>{esc(str(r.size_bp))}</td><td class='num'>{esc(r.split_read_support)}</td>"
                             f"<td>{esc(r.genes_overlapped)}</td></tr>" for r in sv.head(50).itertuples())
                   + "</tbody></table>")
    else:
        sv_html = "<p>No mutant-specific structural variants passed the filters.</p>"

# ------------------------------------------------------------------ plots + peak
PLOT_TEXT = {
    "candidate_map": ("Where are the candidates?",
        "Each grey bar is a chromosome. Orange dots are priority-1 variants, blue ticks priority-2, "
        "grey triangles structural variants that overlap a gene. The shaded box marks the region "
        "with the most homozygous mutant-specific variants."),
    "variant_density": ("Mutant-specific variant density",
        "Number of high-confidence mutant-specific variants in each 1 Mb window (top), and how many of them "
        "are homozygous in the mutant (bottom). A tall, isolated peak means the mutant differs from the "
        "control over a whole region – typically the chromosome segment around a recessive mutation."),
    "homozygosity": ("Homozygosity",
        "Top: fraction of variant sites that are homozygous in each sample. Bottom: mutant minus control. "
        "Red blocks are regions where the mutant has lost variation (both copies identical) but the control "
        "has not – the expected footprint of a recessive mutation."),
    "homozygosity_runs": ("Runs of homozygosity",
        "Runs where both chromosome copies are identical, found with bcftools roh (a hidden Markov model "
        "that weighs every genotype call along the chromosome). Top: share of each 1 Mb window inside a run, "
        "per sample. Bottom: mutant minus control; the shaded block is the strongest stretch where the "
        "mutant is in runs and the control is not. At low coverage runs come out fragmented, so read the "
        "window shares rather than single runs (homozygosity/roh_runs.tsv)."),
    "homozygosity_mapping": ("Homozygosity mapping",
        "At sites where the mutant carries a variant and the control shows the reference letter, this is the "
        "share of mutant reads that still show the reference, i.e. come from another version of the "
        "chromosome. Next to a recessive mutation every mutant chromosome carries the same version, so the "
        "share drops to about 0 (sequencing errors give 1–2%); two different copies give the dashed "
        "genome-typical level. Values in between mean that only some of the mutant DNA carries another "
        "version, which happens when the sample is a pool of fish and some of them carry a crossover. "
        "Top: whole genome. Bottom: the chromosome with the lowest point. The dark band (best interval) "
        "is where the mutation most likely lies; per-window values: homozygosity/mapping_windows.tsv."),
    "manhattan": ("Manhattan plot",
        "Every dot is one variant; the higher it is, the more strongly its reads favour the mutant over the "
        "control (one-sided Fisher exact test). Look for towers of dots rather than single high points. "
        "At low coverage the maximum height is limited, so no genome-wide significance line is drawn."),
    "coverage": ("Read depth and copy number",
        "Top: read depth in 500 kb windows relative to the genome median for each sample. Bottom: "
        "log2 ratio mutant/control. A step down to about −1 suggests a large deletion of one copy in the "
        "mutant; up to about +0.58 a duplication. Dips seen in both samples are repetitive regions."),
    "summary": ("Summary",
        "(a) variant categories, (b) the kinds of protein changes among high-confidence mutant-specific "
        "variants, (c) allele fractions in the mutant, (d) sizes of mutant-specific structural variants."),
}
pkc = pd.read_csv(inp.peak_candidates, sep="\t", dtype=str).fillna("")
MAP_PILL = {"best interval": '<span class="pill p1">best interval</span>',
            "linked": '<span class="pill p2">linked</span>', "outside": "outside"}
if "mapping" in pkc.columns:
    order = {"best interval": 0, "linked": 1, "outside": 2}
    pkc = pkc.assign(_o=pkc["mapping"].map(order), _p=pd.to_numeric(pkc["priority"]),
                     _x=pd.to_numeric(pkc["pos"])).sort_values(["_o", "_p", "_x"]).drop(columns=["_o", "_p", "_x"])
peak_table = ""
if len(pkc):
    peak_table = (
        f'<div class="card fig"><h3>Protein-changing variants inside the peak ({len(pkc)})</h3>'
        '<p class="sub">Every mutant-specific HIGH/MODERATE variant in the peak region, all priorities. Inside a '
        'linked region most variants are homozygous "passengers", so expect several; the causal change is '
        'often the one that best fits the phenotype (check ZFIN). Full list: '
        '<code>mutant_specific/peak_candidates.tsv</code>. '
        + ('<b>Mapping</b>: whether the variant lies in the best interval or the wider linked region of the '
           'homozygosity-mapping plot; the table is sorted with the best interval first.</p>'
           if "mapping" in pkc.columns else "</p>")
        + '<table><thead><tr><th class="num">Priority</th><th>Location</th><th>Gene</th><th>Effect</th>'
        '<th>Change</th><th>Confidence</th>' + ('<th>Mapping</th>' if "mapping" in pkc.columns else '')
        + '</tr></thead><tbody>'
        + "".join(f"<tr><td class='num'>{prio_badge(r.priority)}</td><td>{loc_link(r.chrom, r.pos)}</td>"
                  f"<td>{gene_link(r.gene, r.gene_id)}</td><td>{esc(r.effect)}</td>"
                  f"<td><code>{esc(r.hgvs_p)}</code></td>"
                  f"<td>{'high' if r.category in HIGH_CONF else 'low'}</td>"
                  + (f"<td>{MAP_PILL.get(r.mapping, esc(r.mapping))}</td>" if "mapping" in pkc.columns else "")
                  + "</tr>"
                  for r in pkc.head(150).itertuples())
        + "</tbody></table></div>")

plot_by_name = {p_.rsplit("/", 1)[-1][:-4]: p_ for p_ in inp.plots}
plots_html = ""
for name, (title, text) in PLOT_TEXT.items():
    if name in plot_by_name:
        b64 = base64.b64encode(open(plot_by_name[name], "rb").read()).decode()
        plots_html += (f'<div class="card fig"><h3>{title}</h3><p class="sub">{text}</p>'
                       f'<img alt="{title}" src="data:image/png;base64,{b64}">'
                       f'<div class="sub">plots/{name}.png · plots/{name}.pdf</div></div>')
        if name == [n for n in ("homozygosity", "homozygosity_runs", "homozygosity_mapping") if n in plot_by_name][-1]:
            plots_html += peak_table

peak = json.load(open(inp.peak))
peak_html = ""
roh_peak = peak.get("runs_of_homozygosity")
roh_sentence = ""
if roh_peak:
    agrees = (peak.get("found") and str(roh_peak["chrom"]) == str(peak["chrom"])
              and roh_peak["start_mb"] < peak["end_mb"] and roh_peak["end_mb"] > peak["start_mb"])
    roh_sentence = (
        f' Runs of homozygosity {"agree" if agrees else "point elsewhere"}: the strongest stretch where the '
        f'mutant is homozygous and the control is not is chromosome {esc(str(roh_peak["chrom"]))}, '
        f'{roh_peak["start_mb"]}–{roh_peak["end_mb"]} Mb (mutant in runs over '
        f'{100 * roh_peak["mutant_share"]:.0f}% of it, control {100 * roh_peak["control_share"]:.0f}%).')
mapping = peak.get("homozygosity_mapping")
map_sentence = ""
if mapping:
    mc = mapping.get("candidates", [])
    top_ = [c for c in mc if str(c["priority"]) in ("1", "2")]
    genes_ = list(dict.fromkeys(str(c["label"]) for c in top_))
    map_sentence = (
        f' <b>Homozygosity mapping</b> narrows it to chromosome {esc(str(mapping["chrom"]))}, '
        f'{mapping["start"] / 1e6:.2f}–{mapping["end"] / 1e6:.2f} Mb: there, only '
        f'{100 * mapping["lowest_share"]:.1f}% of mutant reads come from another version of the chromosome '
        f'(genome typical {100 * mapping["typical_share"]:.0f}%; next-best chromosome '
        f'{esc(str(mapping["next_chrom"]))} at {100 * mapping["next_share"]:.0f}%). '
        + (f'It holds {len(mc)} protein-changing mutant-specific variant(s), {len(top_)} of priority 1–2 in '
           f'{", ".join("<b>" + esc(g) + "</b>" for g in genes_)}.' if top_ else
           f'It holds {len(mc)} protein-changing mutant-specific variant(s), none of priority 1–2.'))
if peak.get("found"):
    share = 100 * peak["homozygous_variants_in_peak"] / max(peak["homozygous_variants_total"], 1)
    cands = peak["candidates"]
    clist = ", ".join(
        f"<b>{esc(str(c['label']))}</b> (priority {c['priority']}, {esc(str(c['effect']).split(',')[0])})"
        for c in cands[:12]) + (f" and {len(cands) - 12} more" if len(cands) > 12 else "")
    peak_html = (
        f'<div class="note peak"><b>Strongest genomic signal: chromosome {esc(str(peak["chrom"]))}, '
        f'{peak["start_mb"]}–{peak["end_mb"]} Mb.</b> This region holds {share:.0f}% of all homozygous '
        f'mutant-specific variants (up to {peak["peak_per_mb"]:,} per Mb vs a genome average of '
        f'{peak["genome_mean_per_mb"]}). Regions like this usually surround a recessive mutation, but they can '
        f'also be a stretch the two samples simply inherited differently. '
        + (f'{len(cands)} priority 1–2 variant(s) fall inside it: {clist}. Check these first.'
           if cands else 'No priority 1–2 variants fall inside it; look at priority 3 variants and '
                         'structural variants in this region, or in genes poorly covered by the reads.')
        + ' Priority (1–5) is scored one variant at a time and ignores location, so the top-priority '
          'variants elsewhere in the genome can rank above variants in this region; weigh both. '
          'See the plots and the table of variants inside the peak in section 4.'
        + roh_sentence + map_sentence + '</div>')

# ------------------------------------------------------------------ warnings
warnings = []
if any(not p for p in paired.values()):
    warnings.append("Some samples are <b>single-end</b> (only read 1 or read 2). Alignment and small-variant "
                    "calling work well, but structural-variant detection is much less sensitive without read pairs.")
low = [r["sample"] for r in qc if r["mean_depth"] < 10]
if low:
    warnings.append(f"Average coverage is below 10× for <b>{', '.join(map(esc, low))}</b>. At low coverage one "
                    "copy of a heterozygous variant is often missed, so some 'mutant-specific' variants may in fact "
                    "also be present in the control. Prioritise <i>homozygous</i>, high-confidence hits and "
                    "confirm candidates by PCR/Sanger sequencing.")
if len(summary["controls"]) < 2 or len(mutants) < 2:
    warnings.append("With one sample per group, many mutant-specific variants are ordinary inherited differences "
                    "between the fish sampled, not the causal mutation, unless each sample is a large pool. "
                    "More controls (or pooled siblings) remove most of this background.")
if sibling:
    warnings.append("Comparison mode <b>sibling</b>: candidates must be carried in (nearly) all mutant reads "
                    f"(≥ {cfg['comparison'].get('min_alt_fraction_mutant', 0.9):g}), while the sibling control may "
                    f"carry them in up to {cfg['comparison'].get('max_alt_fraction_control', 0.6):g} of its reads.")
elif summary["controls"]:
    warnings.append("Comparison mode <b>absent</b>: a variant only counts as mutant-specific when the control shows "
                    "no trace of it. If your control is a <i>sibling</i> (or sibling pool) of a recessive mutant, it "
                    "may carry one copy of the mutation and hide it – set <code>mode: sibling</code> in "
                    "config/config.yaml.")
warn_html = "".join(f'<div class="note"><span class="pill warn">! Note</span> {w}</div>' for w in warnings)

n_hc = sum(cats.get(k, 0) for k in HIGH_CONF)
n_genes = len(genes)
pc = {str(k): v for k, v in summary.get("priority_counts", {}).items()}
qc_ok = all(light(r["mean_depth"], 15, 8)[0] != "bad" and light(r["mapped_pct"], 95, 85)[0] != "bad" for r in qc)

tiles = [
    ("Samples", f"{len(names)}", f"{len(summary['controls'])} control · {len(mutants)} mutant"),
    ("Data quality", "Pass" if qc_ok else "Check", "see table below"),
    ("Mutant-enriched variants" if sibling else "Mutant-specific variants", f"{n_hc:,}", f"+ {cats.get('mutant_specific_low_confidence', 0):,} low confidence"),
    ("Top candidates", f"{pc.get('1', 0)} + {pc.get('2', 0)}", "priority 1 + priority 2 variants"),
    ("Candidate genes", f"{n_genes:,}", "with any protein-changing hit"),
]
if n_sv is not None:
    tiles.append(("Structural variants", f"{n_sv:,}", "mutant-specific"))
tile_html = "".join(f'<div class="tile"><div class="tl">{a}</div><div class="tv">{b}</div><div class="sub">{c}</div></div>'
                    for a, b, c in tiles)

CSS = """
:root{--bg:#fcfcfb;--card:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--line:#e4e3df;--accent:#2a78d6;
--series-1:#2a78d6;--series-2:#eb6834;--series-3:#1baf7a;--good:#0a7d33;--good-bg:#e3f4e8;--warn:#9a5b00;--warn-bg:#fdf0d5;
--bad:#b42318;--bad-bg:#fde6e4;--bar:#86b6ef}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#1a1a19;--card:#232322;--ink:#f0efec;--ink2:#c3c2b7;
--line:#3a3936;--accent:#3987e5;--series-1:#3987e5;--series-2:#d95926;--series-3:#199e70;--good:#7fd69a;--good-bg:#173a24;
--warn:#f0c26b;--warn-bg:#3d2f12;--bad:#f59b93;--bad-bg:#43201d;--bar:#256abf}}
:root[data-theme="dark"]{--bg:#1a1a19;--card:#232322;--ink:#f0efec;--ink2:#c3c2b7;--line:#3a3936;--accent:#3987e5;
--series-1:#3987e5;--series-2:#d95926;--series-3:#199e70;--good:#7fd69a;--good-bg:#173a24;--warn:#f0c26b;--warn-bg:#3d2f12;
--bad:#f59b93;--bad-bg:#43201d;--bar:#256abf}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}
main{max-width:1100px;margin:0 auto;padding:24px 16px 64px}h1{font-size:26px;margin:0 0 4px}h2{font-size:19px;margin:36px 0 8px}
a{color:var(--accent)}.sub{color:var(--ink2);font-size:13px}code{font-size:12.5px}
.card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:16px;overflow-x:auto}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin-top:18px}
.tile{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:14px}.tl{color:var(--ink2);font-size:13px}
.tv{font-size:28px;font-weight:650;font-variant-numeric:tabular-nums}
table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:7px 10px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-weight:600;color:var(--ink2);font-size:13px}.num{text-align:right;font-variant-numeric:tabular-nums;white-space:nowrap}
.pill{display:inline-block;font-size:11.5px;padding:1px 7px;border-radius:99px;margin-left:6px;font-weight:600}
.pill.good{color:var(--good);background:var(--good-bg)}.pill.warn{color:var(--warn);background:var(--warn-bg)}
.pill.bad{color:var(--bad);background:var(--bad-bg)}
.pill.p1{color:var(--card);background:var(--accent)}.pill.p2{color:var(--accent);border:1px solid var(--accent)}.pill.na{color:var(--ink2)}
.note{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px 14px;margin:8px 0}
.barcell{width:28%}.bar{height:12px;background:var(--bar);border-radius:0 4px 4px 0;min-width:2px}
svg{width:100%;height:auto;display:block}.grid{stroke:var(--line);stroke-width:1}.tick{fill:var(--ink2);font-size:11px}
.lbl{fill:var(--ink);font-size:12px}.dot{stroke:var(--card);stroke-width:2}.dot:hover{r:6}
.legend{display:flex;gap:16px;font-size:13px;margin-bottom:6px}.key i{display:inline-block;width:12px;height:12px;border-radius:3px;margin-right:6px;vertical-align:-1px}
details summary{cursor:pointer;color:var(--accent)}
.fig{margin-bottom:14px}.fig h3{margin:0 0 4px;font-size:16px}.fig img{width:100%;height:auto;display:block;
background:#fff;border-radius:6px;margin:8px 0 4px}.note.peak{border-left:4px solid var(--series-2)}ul.files li{margin:4px 0}
"""

files = [
    ("mutant_specific/mutant_specific_variants.xlsx", "<b>Start here.</b> Excel workbook of likely damaging mutant-specific variants and candidate genes."),
    ("mutant_specific/candidate_genes.tsv", "One row per gene hit by a protein-changing mutant-specific variant."),
    ("mutant_specific/mutant_specific_variants.tsv.gz", "Every mutant-specific variant, including non-coding ones (compressed text table; often too large for Excel)."),
    ("mutant_specific/mutant_specific_variants.vcf.gz", "Same variants in VCF format (for IGV or other tools)."),
    ("homozygosity/roh_runs.tsv", "Runs of homozygosity in each sample (bcftools roh): start, end, length."),
    ("homozygosity/mapping_windows.tsv", "Homozygosity mapping: share of mutant reads from another chromosome version in each window."),
    ("structural_variants/mutant_specific_sv.tsv", "Large deletions, duplications, inversions found only in the mutant."),
    ("variants/all_samples.annotated.vcf.gz", "All variants in all samples, with gene annotation."),
    ("alignments/", "Aligned reads (BAM). Open in IGV to look at a variant by eye."),
    ("qc/multiqc_report.html", "Detailed quality-control report."),
]
file_html = "".join(f'<li><a href="{f}">{f}</a> – {d}</li>' for f, d in files)

vc = cfg["variant_calling"]
cc = cfg["comparison"]
if sibling:
    methods_compare = (
        f"Mutant and sibling control were compared in sibling mode: a candidate had ≥ {cc['min_alt_reads_mutant']} "
        f"supporting reads in the mutant and an alternate-allele fraction ≥ {cc.get('min_alt_fraction_mutant', 0.9):g}, "
        f"and each control (depth ≥ {cc['min_depth_control']}) had an alternate-allele fraction ≤ "
        f"{cc.get('max_alt_fraction_control', 0.6):g} (mutant-specific if zero, mutant-enriched otherwise).")
else:
    methods_compare = (
        f"A variant was called mutant-specific if the mutant carried ≥ {cc['min_alt_reads_mutant']} supporting reads "
        f"and an alternate genotype, and each control had zero supporting reads with depth ≥ {cc['min_depth_control']}.")
if summary.get("known_variants_checked"):
    methods_compare += (" Variants were matched (same position and allele, after normalisation) to known zebrafish "
                        "variation from Ensembl (dbSNP/EVA).")
if summary.get("protein_check"):
    sp = ", ".join(cfg.get("protein_check", {}).get("species", []))
    methods_compare += (" For missense variants, the zebrafish protein of the annotated transcript was aligned "
                        "(Biopython, BLOSUM62) to the NCBI RefSeq proteins of the same gene and to the Ensembl "
                        f"Compara ortholog in each of {esc(sp)}, and the aligned residues were compared with the "
                        "reference and mutant amino acids.")
methods_roh = ""
if roh_peak is not None or "homozygosity_runs" in plot_by_name:
    rc = cfg.get("homozygosity_runs", {})
    methods_roh = (
        f" Runs of homozygosity were detected with bcftools roh on SNP genotypes (genotype-only mode, "
        f"-G {rc.get('genotype_error_phred', 30)}, allele frequency {rc.get('allele_frequency', 0.4)}, "
        f"AZ→HW transition {rc.get('az_to_hw', 1e-10):g}, hom-ref calls ignored), after setting calls with fewer "
        f"than {rc.get('min_reads', 3)} reads, or heterozygous calls with fewer than "
        f"{rc.get('min_reads_each_allele', 2)} reads of either allele, to missing; the share of each 1 Mb "
        f"window in runs was compared between mutant and control.")
methods_map = ""
if mapping is not None or "homozygosity_mapping" in plot_by_name:
    hm = cfg.get("homozygosity_mapping", {})
    methods_map = (
        f" For homozygosity mapping, SNVs where the mutant had ≥ 2 alternate reads (depth ≥ 4) and the control "
        f"≥ 1 reference read (depth ≥ 3) were grouped in {hm.get('window_kb', 250)} kb windows (≥ "
        f"{hm.get('min_sites', 30)} sites); the share of mutant reads showing the reference allele was computed "
        f"per window and smoothed over 5 windows. The best interval is the run of windows, around the lowest "
        f"smoothed point in the genome, whose share is within {hm.get('core_margin', 0.02):g} of the lowest window.")
page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>WGS Variant Report</title>
<style>{CSS}</style></head><body><main>
<h1>Whole-genome sequencing report</h1>
<div class="sub">Generated {datetime.datetime.now():%Y-%m-%d %H:%M} · reference {esc(cfg['reference']['name'])} (zebrafish)</div>
<div class="tiles">{tile_html}</div>
{peak_html}
{warn_html}

<h2>1. Is the data good enough?</h2>
<p class="sub"><b>Coverage</b> = how many reads cover each DNA letter on average (≥15× ideal, ≥8× usable).
<b>Mapped</b> = reads that matched the zebrafish genome. <b>Duplicates</b> = copies of the same DNA fragment (lower is better).
<b>Genome covered</b> = percent of the genome seen by at least one read.</p>
<div class="card"><table><thead><tr><th>Sample</th><th>Group</th><th>Reads</th><th class="num">Raw reads</th>
<th class="num">Reads kept</th><th class="num">High-quality bases (Q30)</th><th class="num">Mapped</th>
<th class="num">Duplicates</th><th class="num">Mean coverage</th><th class="num">Genome covered</th></tr></thead>
<tbody>{''.join(qc_rows)}</tbody></table></div>

<h2>2. Coverage along the genome</h2>
<p class="sub">Mean coverage on chromosomes 1–25. Lines should be flat; a sample-specific dip can mean a large deletion
or a chromosome loss. Hover a point for its value.</p>
<div class="card">{coverage_svg()}</div>

<h2>3. How variants were classified</h2>
<p class="sub">A <b>variant</b> is a place where the DNA differs from the zebrafish reference genome.
<b>HIGH</b> impact = likely breaks the protein (stop codon, frameshift, splice site);
<b>MODERATE</b> = changes one amino acid (missense).</p>
<div class="card"><table><thead><tr><th>Category</th><th class="num">Variants</th><th></th>
<th class="num">HIGH</th><th class="num">MODERATE</th></tr></thead><tbody>{cat_rows}</tbody></table></div>
{checks_html}

<h2>4. Genome-wide plots</h2>
<p class="sub">Each figure is also saved as PNG (slides) and PDF (papers) in <code>plots/</code>.</p>
{plots_html}

<h2>5. Candidate genes (top 50)</h2>
<p class="sub">Genes with protein-changing variants found only in the mutant, ranked by their best variant's priority (1 = strongest), then by
the number of high-confidence hits. A recessive mutant phenotype usually needs both copies hit (homozygous).
Full list: <code>mutant_specific/candidate_genes.tsv</code>.</p>
<div class="card"><table><thead><tr><th class="num">Priority</th><th>Gene</th><th>Chr</th><th class="num">HIGH</th><th class="num">MODERATE</th>
<th class="num">High confidence</th><th class="num">Homozygous</th><th>Effects</th></tr></thead>
<tbody>{gene_rows or '<tr><td colspan=8>None found.</td></tr>'}</tbody></table></div>

<h2>6. Top-priority variants (priority 1 and 2, up to 100)</h2>
<p class="sub"><b>Priority 1</b>: HIGH impact, high confidence, homozygous in the mutant. <b>Priority 2</b>: HIGH impact high confidence on one copy, or MODERATE impact high confidence homozygous. Priority 3–5 are in the Excel file (5 = probably not the mutation, see its 'warnings' column). <b>Checks</b>: whether the variant is known/novel and how the amino acid compares with other zebrafish and vertebrate proteins. Click a location to see it in the Ensembl genome browser.</p>
<div class="card"><table><thead><tr><th class="num">Priority</th><th>Location</th><th>Gene</th><th>Effect</th><th>Change</th><th>Checks</th>{var_head}</tr></thead>
<tbody>{var_rows or '<tr><td colspan=10>None found.</td></tr>'}</tbody></table></div>

<h2>7. Structural variants (top 50)</h2>
<p class="sub">Large changes (deletions, duplications, inversions, translocations) seen only in the mutant.</p>
<div class="card">{sv_html}</div>

<h2>8. Where are my files?</h2>
<div class="card"><ul class="files">{file_html}</ul></div>

<h2>9. Methods</h2>
<div class="card sub">Reads were trimmed with fastp (adapter and poly-G removal, minimum length {cfg['trimming']['min_length']} bp),
aligned to {esc(cfg['reference']['name'])} with bwa-mem2, and PCR duplicates were marked with samtools markdup.
Small variants were called jointly across all samples with bcftools mpileup/call (MAPQ ≥ {cfg['alignment']['min_mapping_quality']},
base quality ≥ {vc['min_base_quality']}), normalised, filtered to QUAL ≥ {vc['min_qual']} and annotated with snpEff.
{methods_compare}{methods_roh}{methods_map} Structural variants were called with DELLY.</div>
</main></body></html>"""

open(snakemake.output[0], "w").write(page)
