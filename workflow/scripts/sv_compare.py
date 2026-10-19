"""Keep DELLY structural variants found in mutant(s) but not control(s),
and list the genes they overlap."""

import bisect
import subprocess
import sys
from collections import defaultdict

import pandas as pd

sys.stderr = open(snakemake.log[0], "w")
bcf = snakemake.input.bcf
controls = list(snakemake.params.controls)
mutants = list(snakemake.params.mutants)
gtf = snakemake.config["reference"]["dir"] + "/genes.gtf"

# gene intervals per chromosome
genes = defaultdict(list)
with open(gtf) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if f[2] != "gene":
            continue
        attrs = dict(a.strip().split(" ", 1) for a in f[8].strip().rstrip(";").split(";") if a.strip())
        name = attrs.get("gene_name", attrs.get("gene_id", "")).strip('"')
        genes[f[0]].append((int(f[3]), int(f[4]), name))
for c in genes:
    genes[c].sort()
starts = {c: [g[0] for g in v] for c, v in genes.items()}


def overlapping(chrom, start, end, max_list=30):
    if chrom not in genes:
        return []
    hits = []
    i = bisect.bisect_right(starts[chrom], end)
    for gs, ge, name in genes[chrom][max(0, i - 2000):i]:
        if ge >= start and gs <= end:
            hits.append(name)
    return hits[:max_list] + (["..."] if len(hits) > max_list else [])


samples = mutants + controls
fmt = "%CHROM\t%POS\t%ID\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%FILTER\t%INFO/PE\t%INFO/SR[\t%GT]\n"
q = subprocess.run(["bcftools", "query", "-s", ",".join(samples), "-f", fmt, bcf],
                   capture_output=True, text=True, check=True)

rows, keep_ids = [], []
for line in q.stdout.splitlines():
    f = line.split("\t")
    chrom, pos, sv_id, svtype, end, svlen, flt, pe, sr = f[:9]
    gts = dict(zip(samples, f[9:]))
    mut_alt = [g not in ("0/0", "./.", ".") for s, g in gts.items() if s in mutants]
    ctrl_ref = [g == "0/0" for s, g in gts.items() if s in controls]
    if flt != "PASS" or not any(mut_alt) or not all(ctrl_ref):
        continue
    start, end_i = int(pos), int(end) if end != "." else int(pos)
    hit = overlapping(chrom, start, end_i)
    rows.append({
        "chrom": chrom, "start": start, "end": end_i, "type": svtype,
        "size_bp": end_i - start if svtype != "BND" else None,
        "paired_end_support": pe, "split_read_support": sr,
        **{f"{s}_genotype": gts[s] for s in samples},
        "genes_overlapped": ", ".join(hit), "id": sv_id,
    })
    keep_ids.append(sv_id)

df = pd.DataFrame(rows)
df.to_csv(snakemake.output.tsv, sep="\t", index=False)

ids_file = snakemake.output.tsv + ".ids"
with open(ids_file, "w") as fh:
    fh.write("\n".join(keep_ids) + "\n")
subprocess.run(["bcftools", "view", "-i", f"ID=@{ids_file}", "-Oz", "-o",
                snakemake.output.vcf, bcf], check=True)
subprocess.run(["rm", "-f", ids_file])
print(f"{len(rows)} mutant-specific SVs", file=sys.stderr)
