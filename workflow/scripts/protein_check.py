"""Check missense alleles against independent zebrafish proteins and vertebrate orthologs.

Two common false leads in mutant screens:
  * the reference genome (GRCz11) carries a rare allele or an error, and the "mutant"
    allele is simply the normal zebrafish residue -> it is found in RefSeq proteins,
    which are built from zebrafish cDNAs;
  * the mutant residue is the one most other vertebrates have -> unlikely to be damaging.

For each missense variant the zebrafish protein of the annotated transcript is aligned
to (a) the RefSeq proteins of the same gene and (b) the Ensembl Compara ortholog of each
chosen species, and the residue aligned to the variant position is read off.

Used by compare_samples.py (imported, not a Snakemake rule).
"""

import gzip
import re
from collections import defaultdict
from multiprocessing import Pool

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Sec": "U"}
SHORT_NAMES = {"homo_sapiens": "human", "mus_musculus": "mouse", "gallus_gallus": "chicken",
               "xenopus_tropicalis": "frog", "oryzias_latipes": "medaka",
               "takifugu_rubripes": "fugu", "lepisosteus_oculatus": "gar",
               "gasterosteus_aculeatus": "stickleback", "rattus_norvegicus": "rat"}
MISSENSE = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
NOT_BLOSUM = re.compile(r"[^ARNDCQEGHILKMFPSTWYVBZX]")  # e.g. U (selenocysteine), *, O
FLANK = 10              # residues on each side used to judge the local alignment
MIN_LOCAL_IDENTITY = 0.3
MAX_REFSEQ_PER_GENE = 4

# Explanations used in the Excel README and the report.
CHECK_LABELS = {
    "matches_zebrafish_refseq": "mutant amino acid is found in a published zebrafish protein (RefSeq, from cDNA): "
                                "probably a reference-genome error or common allele, not a new mutation",
    "matches_vertebrate_consensus": "mutant amino acid is the one most other vertebrates have: probably harmless",
    "alt_seen_in_orthologs": "mutant amino acid occurs in at least one other vertebrate: may be tolerated",
    "conserved": "reference amino acid is the same in all compared vertebrates: change more likely damaging",
    "variable": "position differs between vertebrates",
    "no_orthologs": "no ortholog could be compared; the mutant amino acid is not in any zebrafish RefSeq protein",
    "not_checked": "not a single amino-acid change, or no ortholog/RefSeq protein could be aligned",
}
DEMOTING_CHECKS = {"matches_zebrafish_refseq", "matches_vertebrate_consensus"}


def parse_missense(hgvs_p):
    """'p.Pro393Ala' -> ('P', 393, 'A'); None for anything other than one substitution."""
    m = MISSENSE.match(hgvs_p or "")
    if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
        return None
    return AA3[m.group(1)], int(m.group(2)), AA3[m.group(3)]


def read_fasta(path, keep):
    seqs, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in keep:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            else:
                buf.append(line.strip())
    if name in keep:
        seqs[name] = "".join(buf)
    return seqs


# ---------------------------------------------------------------- alignment
_ALIGNER = None


def _aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio import Align
        from Bio.Align import substitution_matrices
        _ALIGNER = Align.PairwiseAligner(
            mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
            open_gap_score=-10, extend_gap_score=-0.5)
        _ALIGNER.end_gap_score = 0  # partial proteins and isoforms may overhang
    return _ALIGNER


def aligned_residues(query, target, positions):
    """Residue of `target` aligned to each 0-based query position (None = gap/poor alignment)."""
    query, target = NOT_BLOSUM.sub("X", query.upper()), NOT_BLOSUM.sub("X", target.upper())
    try:
        blocks = _aligner().align(query, target)[0].aligned
    except (ValueError, OverflowError):  # unusual sequence: treat as not aligned
        return dict.fromkeys(positions)
    q2t = {}
    for (qs, qe), (ts, _) in zip(blocks[0], blocks[1]):
        for q in range(qs, qe):
            q2t[q] = ts + q - qs
    out = {}
    for p in positions:
        if p not in q2t:
            out[p] = None
            continue
        same = sum(1 for k in range(p - FLANK, p + FLANK + 1)
                   if k != p and k in q2t and query[k] == target[q2t[k]])
        out[p] = target[q2t[p]] if same / (2 * FLANK) >= MIN_LOCAL_IDENTITY else None
    return out


# ---------------------------------------------------------------- per-transcript work
_DATA = {}


def _check_transcript(job):
    transcript, gene_id, variants = job
    zf = _DATA["zf"].get(transcript)
    results = {}
    usable = []
    for key, (ref, pos, alt) in variants:
        if zf and pos <= len(zf) and zf[pos - 1] == ref:
            usable.append((key, ref, pos - 1, alt))
        else:
            results[key] = {"allele_check": "not_checked", "refseq_residues": "",
                            "ortholog_residues": ""}
    if not usable:
        return results
    positions = [p for _, _, p, _ in usable]

    refseq = {}
    for acc in _DATA["refseq_by_gene"].get(gene_id, []):
        if acc in _DATA["refseq_pep"]:
            refseq[acc] = aligned_residues(zf, _DATA["refseq_pep"][acc], positions)
    orth = {}
    for sp, prot in _DATA["orth_by_gene"].get(gene_id, []):
        if prot in _DATA["orth_pep"]:
            orth[sp] = aligned_residues(zf, _DATA["orth_pep"][prot], positions)

    for key, ref, p, alt in usable:
        rs = {acc: res[p] for acc, res in refseq.items() if res[p]}
        os_ = {sp: res[p] for sp, res in orth.items() if res[p]}
        n, n_ref = len(os_), sum(r == ref for r in os_.values())
        n_alt = sum(r == alt for r in os_.values())
        if alt in rs.values():
            check = "matches_zebrafish_refseq"
        elif n >= 3 and n_alt > n_ref and 2 * n_alt >= n:
            check = "matches_vertebrate_consensus"
        elif n_alt:
            check = "alt_seen_in_orthologs"
        elif n >= 3 and n_ref == n:
            check = "conserved"
        elif n:
            check = "variable"
        elif rs:
            check = "no_orthologs"
        else:
            check = "not_checked"
        results[key] = {
            "allele_check": check,
            "refseq_residues": ", ".join(f"{a}:{r}" for a, r in rs.items()),
            "ortholog_residues": ", ".join(f"{SHORT_NAMES.get(sp, sp)}:{r}" for sp, r in os_.items()),
        }
    return results


def run_checks(variants, refs, species, threads=1):
    """variants: iterable of (key, transcript, gene_id, hgvs_p).

    Returns {key: {"allele_check", "refseq_residues", "ortholog_residues"}} for missense variants.
    """
    jobs = defaultdict(list)
    genes = {}
    for key, transcript, gene_id, hgvs_p in variants:
        mis = parse_missense(hgvs_p)
        if mis:
            transcript = transcript.split(".")[0]
            jobs[transcript].append((key, mis))
            genes[transcript] = gene_id
    if not jobs:
        return {}

    want_genes = set(genes.values())
    refseq_by_gene = defaultdict(list)
    with open(refs["refseq"]) as fh:
        next(fh)
        rows = [line.rstrip("\n").split("\t") for line in fh]
    # prefer the RefSeq protein linked to the same transcript, then curated (NP_) ones
    rows.sort(key=lambda r: (r[1] not in jobs, not r[2].startswith("NP_")))
    for gene_id, _, acc, _ in rows:
        if gene_id in want_genes and acc not in refseq_by_gene[gene_id] \
                and len(refseq_by_gene[gene_id]) < MAX_REFSEQ_PER_GENE:
            refseq_by_gene[gene_id].append(acc)
    orth_by_gene = defaultdict(list)
    with open(refs["orthologs"]) as fh:
        next(fh)
        for line in fh:
            gene_id, sp, prot = line.split("\t")[:3]
            if gene_id in want_genes and sp in species:
                orth_by_gene[gene_id].append((sp, prot.split(".")[0]))

    _DATA.update(
        zf=read_fasta(refs["zf_pep"], set(jobs)),
        refseq_by_gene=refseq_by_gene,
        refseq_pep=read_fasta(refs["refseq_pep"], {a for v in refseq_by_gene.values() for a in v}),
        orth_by_gene=orth_by_gene,
        orth_pep=read_fasta(refs["orth_pep"], {p for v in orth_by_gene.values() for _, p in v}),
    )
    work = [(t, genes[t], v) for t, v in jobs.items()]
    results = {}
    if threads > 1:
        with Pool(threads) as pool:  # forked workers share _DATA
            for res in pool.imap_unordered(_check_transcript, work, chunksize=20):
                results.update(res)
    else:
        for job in work:
            results.update(_check_transcript(job))
    return results
