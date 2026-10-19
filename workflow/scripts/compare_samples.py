"""Classify annotated variants as mutant-specific, shared or control-only.

Input : jointly called, snpEff-annotated, biallelic VCF (all samples); the ID column
        holds the rsID of known Ensembl/dbSNP variants when that annotation is on.
Output: classified tables, mutant-specific VCF, candidate gene list,
        an Excel workbook for non-programmers and a JSON summary.

Two comparison modes (config comparison: mode):
  absent  - control is unrelated wild type: the variant must be absent from every control.
  sibling - control is a sibling or sibling pool of a recessive mutant: the mutant must
            carry the variant in (nearly) all reads; the control may carry it in up to
            max_alt_fraction_control of its reads (carrier siblings).
"""

import gzip
import json
import subprocess
import sys
from collections import Counter, defaultdict

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
sys.path.insert(0, snakemake.scriptdir)
import protein_check  # noqa: E402

vcf_in = snakemake.input.vcf
out = snakemake.output
par = snakemake.params
controls = list(par.controls)
mutants = list(par.mutants)
min_alt = int(par.min_alt)
min_ctrl_depth = int(par.min_ctrl_depth)
min_depth = int(par.min_depth)
mutant_rule = str(par.rule).lower()
mode = str(par.mode).lower()
min_af_mut = float(par.min_af_mutant)
max_af_ctrl = float(par.max_af_control)
demote_flagged = bool(par.demote_flagged)
known_on = bool(par.known)
protein_refs = dict(par.protein_refs)

if not mutants:
    sys.exit("No sample has group 'mutant' in the sample sheet.")

# categories written to the mutant-specific outputs; the first two are "high confidence"
HIGH_CONF = ["mutant_specific", "mutant_enriched"]
CANDIDATE_CATS = HIGH_CONF + ["mutant_specific_low_confidence"]
ALL_CATS = CANDIDATE_CATS + ["mutant_heterozygous", "shared", "control_only", "low_evidence"]
IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3, "": 4}
CHECK_RANK = {"conserved": 0, "variable": 1, "no_orthologs": 2, "not_checked": 2, "alt_seen_in_orthologs": 3,
              "matches_vertebrate_consensus": 4, "matches_zebrafish_refseq": 5}
ANN_FIELDS = ["effect", "impact", "gene", "gene_id", "transcript", "biotype",
              "hgvs_c", "hgvs_p"]


def parse_ann(info):
    """Return the most severe snpEff annotation (snpEff lists it first)."""
    for item in info.split(";"):
        if item.startswith("ANN="):
            first = item[4:].split(",")[0].split("|")
            first += [""] * (16 - len(first))
            return {
                "effect": first[1].replace("&", ", "),
                "impact": first[2],
                "gene": first[3],
                "gene_id": first[4],
                "transcript": first[6],
                "biotype": first[7],
                "hgvs_c": first[9],
                "hgvs_p": first[10],
            }
    return dict.fromkeys(ANN_FIELDS, "")


def sample_values(fmt_keys, field):
    vals = dict(zip(fmt_keys, field.split(":")))
    gt = vals.get("GT", "./.")
    ad = vals.get("AD", ".").split(",")
    ref_n = int(ad[0]) if ad and ad[0].isdigit() else 0
    alt_n = int(ad[1]) if len(ad) > 1 and ad[1].isdigit() else 0
    depth = ref_n + alt_n
    has_alt = any(a not in ("0", ".") for a in gt.replace("|", "/").split("/"))
    if not has_alt:
        zyg = "ref" if gt not in ("./.", ".") else "no call"
    else:
        alleles = gt.replace("|", "/").split("/")
        zyg = "homozygous" if len(set(alleles)) == 1 else "heterozygous"
    return gt, ref_n, alt_n, depth, has_alt, zyg


def combine(flags):
    return all(flags) if mutant_rule == "all" else any(flags)


def classify_absent(mut_has, mut_af, ctrl_alt, ctrl_depth):
    ctrl_has_alt = any(a > 0 for a in ctrl_alt)
    if combine(mut_has):
        if not controls:
            return "mutant_specific"
        if ctrl_has_alt:
            return "shared"
        if all(d >= min_ctrl_depth for d in ctrl_depth):
            return "mutant_specific"
        return "mutant_specific_low_confidence"
    if ctrl_has_alt:
        return "control_only"
    return "low_evidence"


def classify_sibling(mut_has, mut_af, ctrl_alt, ctrl_depth):
    ctrl_has_alt = any(a > 0 for a in ctrl_alt)
    if combine([h and af >= min_af_mut for h, af in zip(mut_has, mut_af)]):
        if not controls:
            return "mutant_specific"
        covered = [d >= min_ctrl_depth for d in ctrl_depth]
        if any(c and a / d > max_af_ctrl for c, a, d in zip(covered, ctrl_alt, ctrl_depth)):
            return "shared"  # a well-covered control looks homozygous: not the recessive mutation
        if not all(covered):
            return "mutant_specific_low_confidence"
        return "mutant_enriched" if ctrl_has_alt else "mutant_specific"
    if any(mut_has):
        return "shared" if ctrl_has_alt else "mutant_heterozygous"
    if ctrl_has_alt:
        return "control_only"
    return "low_evidence"


classify = classify_sibling if mode == "sibling" else classify_absent


def pooled_fraction(alt, depth):
    return sum(alt) / sum(depth) if sum(depth) else None


view = subprocess.Popen(["bcftools", "view", vcf_in], stdout=subprocess.PIPE, text=True)
vcf_out = subprocess.Popen(["bgzip", "-c"], stdin=subprocess.PIPE,
                           stdout=open(out.vcf, "wb"), text=True)

cols = (["category", "chrom", "pos", "ref", "alt", "type", "qual", "known_id", "gene", "gene_id",
         "effect", "impact", "hgvs_c", "hgvs_p", "transcript", "biotype", "af_difference"]
        + [c for s in mutants + controls for c in
           (f"{s}_genotype", f"{s}_ref_reads", f"{s}_alt_reads", f"{s}_alt_fraction")])

all_out = gzip.open(out.all_tsv, "wt", compresslevel=3)
all_out.write("\t".join(cols) + "\n")
sample_idx = {}
counts = Counter()
impact_counts = defaultdict(Counter)
mut_rows = []

for line in view.stdout:
    if line.startswith("##"):
        vcf_out.stdin.write(line)
        continue
    if line.startswith("#CHROM"):
        vcf_out.stdin.write(line)
        header = line.rstrip("\n").split("\t")
        sample_idx = {s: header.index(s) for s in controls + mutants}
        continue
    f = line.rstrip("\n").split("\t")
    fmt = f[8].split(":")
    ann = parse_ann(f[7])
    rec = {
        "chrom": f[0], "pos": int(f[1]), "ref": f[3], "alt": f[4],
        "type": "SNV" if len(f[3]) == len(f[4]) == 1 else "indel",
        "qual": round(float(f[5]), 1) if f[5] != "." else None,
        "known_id": "" if f[2] == "." else f[2],
        **ann,
    }
    mut_has, mut_af, mut_alt, mut_depth, ctrl_alt, ctrl_depth = [], [], [], [], [], []
    for s in mutants + controls:
        gt, ref_n, alt_n, depth, has_alt, zyg = sample_values(fmt, f[sample_idx[s]])
        rec[f"{s}_genotype"] = zyg
        rec[f"{s}_ref_reads"] = ref_n
        rec[f"{s}_alt_reads"] = alt_n
        rec[f"{s}_alt_fraction"] = round(alt_n / depth, 2) if depth else None
        if s in mutants:
            mut_has.append(has_alt and alt_n >= min_alt and depth >= min_depth)
            mut_af.append(alt_n / depth if depth else 0.0)
            mut_alt.append(alt_n)
            mut_depth.append(depth)
        else:
            ctrl_alt.append(alt_n)
            ctrl_depth.append(depth)
    m_frac = pooled_fraction(mut_alt, mut_depth)
    c_frac = pooled_fraction(ctrl_alt, ctrl_depth) if controls else 0.0
    rec["af_difference"] = round(m_frac - c_frac, 2) if m_frac is not None and c_frac is not None else None
    cat = classify(mut_has, mut_af, ctrl_alt, ctrl_depth)
    rec["category"] = cat
    counts[cat] += 1
    impact_counts[cat][ann["impact"] or "NONE"] += 1
    all_out.write("\t".join("" if rec[c] is None else str(rec[c]) for c in cols) + "\n")
    if cat in CANDIDATE_CATS:
        vcf_out.stdin.write(line)
        mut_rows.append(rec)

view.wait()
all_out.close()
vcf_out.stdin.close()
vcf_out.wait()
if view.returncode != 0:
    sys.exit("bcftools view failed")
subprocess.run(["tabix", "-f", "-p", "vcf", out.vcf], check=True)

mut_df = pd.DataFrame(mut_rows, columns=cols)

# --- protein check of amino-acid changes (reference error / vertebrate consensus)
check_cols = ["allele_check", "refseq_residues", "ortholog_residues"]
for c in check_cols:
    mut_df[c] = ""
if protein_refs:
    coding = mut_df[mut_df["impact"].isin(["HIGH", "MODERATE"])]
    checks = protein_check.run_checks(
        zip(coding.index, coding["transcript"], coding["gene_id"], coding["hgvs_p"]),
        protein_refs, list(par.species), threads=snakemake.threads)
    for c in check_cols:
        mut_df[c] = mut_df.index.map(lambda i: checks.get(i, {}).get(c, ""))
    print(f"protein check: {len(checks):,} missense variants", file=sys.stderr)


def warnings(r):
    w = []
    if r["known_id"]:
        w.append(f"known variant ({r['known_id']})")
    if r["allele_check"] == "matches_zebrafish_refseq":
        w.append("mutant amino acid is in a published zebrafish protein")
    elif r["allele_check"] == "matches_vertebrate_consensus":
        w.append("mutant amino acid is the usual vertebrate one")
    return "; ".join(w)


mut_df["warnings"] = mut_df.apply(warnings, axis=1) if len(mut_df) else ""
mut_df["novel"] = mut_df["known_id"].eq("") if known_on else ""


def priority(df):
    """1 = strongest candidate ... 5 = weakest; 0 = not protein-changing.

    1: HIGH impact, high confidence, homozygous in a mutant
    2: HIGH impact high confidence (heterozygous), or MODERATE high confidence homozygous
    3: other HIGH/MODERATE high confidence
    4: HIGH/MODERATE low confidence
    5: HIGH/MODERATE but probably not the mutation (known variant, or the mutant amino acid
       is found in a published zebrafish protein / is the usual vertebrate one)
    """
    hc = df["category"].isin(HIGH_CONF)
    hom = df[[f"{s}_genotype" for s in mutants]].eq("homozygous")
    hom = hom.all(axis=1) if mutant_rule == "all" else hom.any(axis=1)
    high, mod = df["impact"] == "HIGH", df["impact"] == "MODERATE"
    p = pd.Series(0, index=df.index)
    p[(high | mod) & ~hc] = 4
    p[(high | mod) & hc] = 3
    p[(high & hc) | (mod & hc & hom)] = 2
    p[high & hc & hom] = 1
    if demote_flagged:
        p[(high | mod) & (df["warnings"] != "")] = 5
    return p


mut_df.insert(0, "priority", priority(mut_df))
mut_df["_p"] = mut_df["priority"].replace(0, 9)
mut_df["_rank"] = mut_df["impact"].map(IMPACT_RANK).fillna(4)
mut_df["_check"] = mut_df["allele_check"].map(CHECK_RANK).fillna(2)
mut_df = (mut_df.sort_values(["_p", "_rank", "_check", "chrom", "pos"])
          .drop(columns=["_p", "_rank", "_check"]))
# keep the new judgement columns next to the protein change
front = ["priority", "category", "chrom", "pos", "ref", "alt", "type", "qual", "gene", "gene_id",
         "effect", "impact", "hgvs_c", "hgvs_p", "warnings", "novel", "known_id", "allele_check",
         "refseq_residues", "ortholog_residues", "af_difference"]
mut_df = mut_df[front + [c for c in mut_df.columns if c not in front]]
mut_df.to_csv(out.mut_tsv, sep="\t", index=False, compression="gzip")

# --- candidate genes: genes hit by HIGH/MODERATE mutant-specific variants
damaging = mut_df[mut_df["impact"].isin(["HIGH", "MODERATE"])].copy()
genes = []
for (gene, gene_id), g in damaging.groupby(["gene", "gene_id"], dropna=False):
    genes.append({
        "best_priority": int(g["priority"].min()),
        "gene": gene,
        "gene_id": gene_id,
        "high_impact_variants": int((g["impact"] == "HIGH").sum()),
        "moderate_impact_variants": int((g["impact"] == "MODERATE").sum()),
        "high_confidence_variants": int(g["category"].isin(HIGH_CONF).sum()),
        "homozygous_in_mutant": int(g[[f"{s}_genotype" for s in mutants]]
                                    .eq("homozygous").any(axis=1).sum()),
        "unflagged_variants": int((g["warnings"] == "").sum()),
        "conserved_missense": int((g["allele_check"] == "conserved").sum()),
        "chrom": g["chrom"].iloc[0],
        "first_pos": int(g["pos"].min()),
        "effects": "; ".join(sorted(set(g["effect"]))),
        "protein_changes": "; ".join(p for p in g["hgvs_p"] if p)[:500],
    })
genes_df = pd.DataFrame(genes, columns=[
    "best_priority", "gene", "gene_id", "high_impact_variants", "moderate_impact_variants",
    "high_confidence_variants", "homozygous_in_mutant", "unflagged_variants", "conserved_missense",
    "chrom", "first_pos", "effects", "protein_changes"])
if not genes_df.empty:
    genes_df = genes_df.sort_values(
        ["best_priority", "high_confidence_variants", "high_impact_variants",
         "homozygous_in_mutant", "moderate_impact_variants"],
        ascending=[True, False, False, False, False])
genes_df.to_csv(out.genes, sep="\t", index=False)

# --- Excel workbook for people who do not use the command line
if mode == "sibling":
    intro = [
        "This workbook lists DNA variants carried in (nearly) all mutant reads, compared with sibling control(s).",
        "Comparison mode 'sibling': siblings may carry one copy, so the control may show the variant in some reads.",
    ]
else:
    intro = ["This workbook lists DNA variants found in the mutant sample(s) but not in the control(s)."]
readme = pd.DataFrame({"Read me": intro + [
    "Sheet 'Likely damaging': variants that change a protein (HIGH or MODERATE impact), most severe first.",
    "Sheet 'Candidate genes': one row per gene hit by a likely damaging variant.",
    "Sheet 'Summary': how many variants fell into each category.",
    "priority 1 = strongest candidate: HIGH impact, high confidence, homozygous in the mutant.",
    "priority 2 = HIGH impact high confidence (one copy), or MODERATE impact high confidence homozygous.",
    "priority 3 = other high-confidence protein-changing variants. priority 4 = low confidence.",
    "priority 5 = probably not the mutation: see the 'warnings' column.",
    "category = mutant_specific: present in mutant, control well covered and shows no trace of it.",
    "category = mutant_enriched (sibling mode): in all mutant reads, in some control reads, as expected for carrier siblings.",
    "category = mutant_specific_low_confidence: present in mutant, but the control had too few reads to be sure.",
    "impact HIGH: e.g. stop gained, frameshift, splice site. MODERATE: e.g. missense (amino-acid change).",
    "warnings: reasons to doubt the variant. A new mutation (e.g. from ENU) is almost never a known variant.",
    "novel = TRUE: not in the Ensembl/dbSNP list of known zebrafish variants. known_id: its rs number if known.",
    "allele_check (amino-acid changes only):",
] + [f"   {k}: {v}" for k, v in protein_check.CHECK_LABELS.items()] + [
    "refseq_residues: amino acid at this position in published zebrafish proteins (NCBI RefSeq).",
    "ortholog_residues: amino acid at this position in the same protein of other vertebrates.",
    "af_difference: share of mutant reads with the variant minus share of control reads (1 = only in mutant).",
    "<sample>_genotype: heterozygous (one copy) or homozygous (both copies) in that sample.",
    "<sample>_alt_reads: number of sequencing reads supporting the variant.",
    "Full list including non-coding variants: mutant_specific_variants.tsv.gz (compressed text table; usually too large for Excel).",
    f"Settings: mode={mode}, min_alt_reads_mutant={min_alt}, min_depth_control={min_ctrl_depth}, "
    f"min_depth={min_depth}, mutant_rule={mutant_rule}"
    + (f", min_alt_fraction_mutant={min_af_mut}, max_alt_fraction_control={max_af_ctrl}"
       if mode == "sibling" else ""),
]})
summary_df = pd.DataFrame(
    [{"category": c, "variants": counts.get(c, 0),
      **{i: impact_counts[c].get(i, 0) for i in ["HIGH", "MODERATE", "LOW", "MODIFIER"]}}
     for c in ALL_CATS if c in counts or c in ("mutant_specific", "shared")])

EXCEL_MAX = 1_000_000
with pd.ExcelWriter(out.xlsx, engine="openpyxl") as xw:
    readme.to_excel(xw, sheet_name="README", index=False)
    damaging.head(EXCEL_MAX).to_excel(xw, sheet_name="Likely damaging", index=False)
    genes_df.head(EXCEL_MAX).to_excel(xw, sheet_name="Candidate genes", index=False)
    summary_df.to_excel(xw, sheet_name="Summary", index=False)
    for ws in xw.book.worksheets:
        ws.freeze_panes = "A2"
        for col in ws.columns:
            width = max(len(str(c.value or "")) for c in list(col)[:200])
            ws.column_dimensions[col[0].column_letter].width = min(max(width + 2, 8), 60)

hc_damaging = damaging[damaging["category"].isin(HIGH_CONF)]
summary = {
    "mode": mode,
    "controls": controls,
    "mutants": mutants,
    "total_variants": int(sum(counts.values())),
    "categories": {k: int(v) for k, v in counts.items()},
    "impact_by_category": {k: dict(v) for k, v in impact_counts.items()},
    "mutant_specific_by_type": mut_df.groupby(["category", "type"]).size()
                               .rename("n").reset_index().to_dict("records"),
    "candidate_genes": int(len(genes_df)),
    "priority_counts": {int(k): int(v) for k, v in mut_df["priority"].value_counts().items() if k},
    "known_variants_checked": known_on,
    "protein_check": bool(protein_refs),
    "protein_changing_high_confidence": {
        "total": int(len(hc_damaging)),
        "known": int((hc_damaging["known_id"] != "").sum()),
        "allele_check": {k: int(v) for k, v in hc_damaging["allele_check"].replace("", "not_missense")
                         .value_counts().items()},
    },
    "settings": {"mode": mode, "min_alt_reads_mutant": min_alt, "min_depth_control": min_ctrl_depth,
                 "min_depth": min_depth, "mutant_rule": mutant_rule,
                 "min_alt_fraction_mutant": min_af_mut, "max_alt_fraction_control": max_af_ctrl,
                 "demote_flagged": demote_flagged},
}
with open(out.summary, "w") as fh:
    json.dump(summary, fh, indent=2)
print("done", counts, file=sys.stderr)
