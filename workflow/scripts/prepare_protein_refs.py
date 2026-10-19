"""Download the protein sequences used to sanity-check missense variants (runs once).

Outputs (all in reference/protein/):
  zebrafish_pep.fa.gz      Ensembl zebrafish proteins, header = transcript ID
  orthologs.tsv            best Ensembl Compara ortholog per zebrafish gene and species
  ortholog_pep.fa.gz       protein sequences of those orthologs, header = protein ID
  refseq.tsv               Ensembl gene -> RefSeq protein accessions (NCBI, cDNA-based)
  refseq_pep.fa.gz         RefSeq zebrafish proteins, header = accession without version
"""

import csv
import gzip
import io
import re
import sys
import urllib.request

sys.stderr = open(snakemake.log[0], "w", buffering=1)  # line-buffered: progress visible in the log

ftp = snakemake.params.ensembl_ftp.rstrip("/")
species = list(snakemake.params.species)
out = snakemake.output


def fetch(url):
    """Open a URL as a text stream, decompressing .gz on the fly."""
    print("downloading", url, file=sys.stderr)
    resp = urllib.request.urlopen(url, timeout=600)
    if url.endswith(".gz"):
        resp = gzip.GzipFile(fileobj=resp)
    return io.TextIOWrapper(resp, encoding="utf-8", errors="replace")


def read_fasta(handle):
    name, seq = None, []
    for line in handle:
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].rstrip("\n"), []
        else:
            seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def pep_url(sp):
    """Ensembl FASTA file names include the assembly, so look it up in the listing."""
    listing = fetch(f"{ftp}/fasta/{sp}/pep/").read()
    m = re.search(r'href="([^"]+\.pep\.all\.fa\.gz)"', listing)
    if not m:
        sys.exit(f"No protein FASTA found for species '{sp}' at {ftp}/fasta/{sp}/pep/")
    return f"{ftp}/fasta/{sp}/pep/{m.group(1)}"


def strip_version(x):
    return x.split(".")[0]


# --- zebrafish proteins, keyed by transcript (snpEff reports the transcript)
with gzip.open(out.zf_pep, "wt") as fh:
    for header, seq in read_fasta(fetch(pep_url("danio_rerio"))):
        m = re.search(r"transcript:(\S+)", header)
        if m:
            fh.write(f">{strip_version(m.group(1))}\n{seq}\n")

# --- orthologs: keep the most similar one per zebrafish gene and species
best = {}
reader = csv.DictReader(fetch(f"{ftp}/tsv/ensembl-compara/homologies/danio_rerio/"
                              f"Compara.{snakemake.params.release}.protein_default.homologies.tsv"),
                        delimiter="\t")
for r in reader:
    sp = r["homology_species"]
    if sp not in species or not r["homology_type"].startswith("ortholog"):
        continue
    key = (r["gene_stable_id"], sp)
    ident = float(r["homology_identity"])
    if key not in best or ident > best[key][1]:
        best[key] = (r["homology_protein_stable_id"], ident, r["homology_type"])
with open(out.orthologs, "w") as fh:
    fh.write("gene_id\tspecies\tprotein_id\tidentity\thomology_type\n")
    for (gene, sp), (prot, ident, htype) in sorted(best.items()):
        fh.write(f"{gene}\t{sp}\t{prot}\t{ident}\t{htype}\n")
print(f"{len(best):,} gene-species ortholog pairs", file=sys.stderr)
missing = [sp for sp in species if not any(k[1] == sp for k in best)]
if missing:
    sys.exit(f"No zebrafish orthologs found for {', '.join(missing)} in Ensembl Compara. "
             "Check the spelling in protein_check: species (Ensembl names such as homo_sapiens), "
             "or choose another species.")

wanted = {strip_version(p) for p, _, _ in best.values()}
with gzip.open(out.orth_pep, "wt") as fh:
    for sp in species:
        for header, seq in read_fasta(fetch(pep_url(sp))):
            pid = strip_version(header.split()[0])
            if pid in wanted:
                fh.write(f">{pid}\n{seq}\n")

# --- RefSeq: proteins derived from zebrafish cDNAs, independent of the GRCz11 sequence
n = 0
with open(out.refseq, "w") as fh:
    fh.write("gene_id\ttranscript_id\trefseq_protein\tsource\n")
    for r in csv.DictReader(fetch(snakemake.params.refseq_xref_url), delimiter="\t"):
        if r["db_name"] in ("RefSeq_peptide", "RefSeq_peptide_predicted"):
            fh.write(f"{r['gene_stable_id']}\t{r['transcript_stable_id']}\t"
                     f"{strip_version(r['xref'])}\t{r['db_name']}\n")
            n += 1
print(f"{n:,} Ensembl-RefSeq protein links", file=sys.stderr)

with gzip.open(out.refseq_pep, "wt") as fh:
    for header, seq in read_fasta(fetch(snakemake.params.refseq_pep_url)):
        fh.write(f">{strip_version(header.split()[0])}\n{seq}\n")
print("done", file=sys.stderr)
