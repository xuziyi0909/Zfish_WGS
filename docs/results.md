---
title: Understanding your results
nav_order: 6
---

# Understanding your results
{: .no_toc }

1. TOC
{:toc}

---

## Start with these two files

| File | Open with | What it is |
|---|---|---|
| `results/report.html` | any web browser | One-page summary: data quality, candidate genes, top variants |
| `results/mutant_specific/mutant_specific_variants.xlsx` | Excel | Your candidate mutations, ready to sort and filter |

## The report, section by section

### 1. Is the data good enough?

Each sample gets a ✓ **Good**, **! Caution** or ✕ **Poor** label for each measure:

| Measure | Good | What it means if it is low |
|---|---|---|
| Reads kept | ≥ 90% | Many reads were too short or low quality after cleaning |
| High-quality bases (Q30) | ≥ 85% | The sequencing run itself was poor |
| Mapped | ≥ 95% | Contamination, or the wrong species |
| Duplicates | ≤ 15% | The library had too little starting DNA |
| Mean coverage | ≥ 15× | Too little sequencing; variants may be missed |
| Genome covered | ≥ 95% | Parts of the genome got no reads at all |

### 2. Coverage along the genome

Each line is one sample. A flat line is normal. If one sample drops on a single
chromosome, that can point to a large deletion or a missing chromosome.

### 3. How variants were classified

Every variant ends up in exactly one group:

| Category | Meaning | Is it a candidate? |
|---|---|---|
| **Mutant-specific (high confidence)** | Seen in the mutant. The control was well covered at that spot and showed no sign of it. | **Yes, focus here** |
| **Mutant-enriched (high confidence)** | *Sibling mode only.* In (nearly) all mutant reads and in some control reads: what a recessive mutation looks like when control siblings carry one copy. | **Yes, focus here** |
| **Mutant-specific (low confidence)** | Seen in the mutant, but the control had too few reads there to be sure it's absent. | Maybe |
| Mutant, not homozygous | *Sibling mode only.* In some but not all mutant reads. | No, for a recessive mutation |
| Shared | Seen in both mutant and control | No: it's background |
| Control only | Seen only in the control | No |
| Low evidence | Too few reads to decide | No |

### 4. Candidate genes

Genes that carry at least one **protein-changing** (HIGH or MODERATE impact)
mutant-specific variant. Every such variant gets a **priority** from 1 (strongest) to 5:

| Priority | Rule | Why |
|---|---|---|
| **1** | HIGH impact + high confidence + **homozygous** in the mutant | Breaks the protein on both copies: the classic recessive mutant |
| **2** | HIGH impact + high confidence on one copy, **or** MODERATE impact + high confidence + homozygous | Strong, but one piece of evidence is weaker |
| **3** | Any other high-confidence protein-changing variant | Worth a look if priorities 1–2 don't explain the phenotype |
| **4** | Protein-changing but **low confidence** (control poorly covered) | Often background; check in IGV first |
| **5** | Protein-changing but **flagged**: a known variant, or an amino-acid change that matches a published zebrafish protein or the usual vertebrate residue (see below) | Probably not the mutation |

Genes are ranked by their best variant's priority, then by the number of
high-confidence hits.

Each gene name links to its **Ensembl** and **ZFIN** pages. Use them to check
whether the gene's known expression or phenotype matches your mutant.

### 5. Top-priority variants

All priority 1 and 2 variants (up to 100).

Click a location to open it in the Ensembl genome browser. The columns for each
sample show its genotype and **supporting reads / total reads**. For example,
`homozygous (6/6)` means all 6 reads at that spot carry the variant.

### Genome-wide plots

All plots are also saved in `results/plots/` as PNG (for slides) and PDF (for papers).

| Plot | What to look for |
|---|---|
| **Candidate map** | Where the priority 1–2 variants and gene-hitting structural variants sit on the chromosomes. The shaded box is the homozygosity peak. |
| **Variant density** | A tall, isolated tower of mutant-specific variants, especially homozygous ones. That's the chromosome segment the mutant inherited differently from the control. For a recessive mutation, the causal change is usually inside it. |
| **Homozygosity** | Red blocks in the difference track are regions where the mutant has lost variation but the control hasn't. That's the footprint of a recessive mutation. |
| **Runs of homozygosity** | Found with `bcftools roh`, a hidden Markov model that weighs every genotype call along the chromosome. The top track is the share of each 1 Mb window that lies in a run, per sample; both samples share many runs from the lab strain's background. The red block in the difference track is where the mutant is in runs and the control isn't. At low coverage single runs are fragmented, so read the window shares, not individual runs. |
| **Manhattan plot** | Each dot is one variant; higher means its reads favour the mutant more strongly. Look for a *tower* of dots rather than single tall points. At low coverage the maximum height is limited. |
| **Coverage / copy number** | A step to about −1 in the bottom track suggests the mutant lost one copy of a large region; about +0.58 suggests a duplication. Dips in *both* samples are just hard-to-sequence regions. |
| **Summary** | Variant categories, kinds of protein changes, allele fractions (0.5 = one copy, 1.0 = both), structural-variant sizes. |

{: .important-title }
> Priority vs. the peak
>
> **Priority** (1–5) is scored one variant at a time: how damaging, how confident,
> homozygous or not. It **ignores location**. The **peak** is regional evidence from
> thousands of variants. A priority-1 variant far from the peak can be a background
> hit, and the causal variant inside the peak may be "only" a missense change (priority 2).
> Weigh both. Variants in the peak are tagged **in peak** in the report and listed in
> `mutant_specific/peak_candidates.tsv`.

### 6. Structural variants

Large deletions (DEL), duplications (DUP), inversions (INV), insertions (INS) and
translocations (BND) found only in the mutant, together with the genes they overlap.

## The Excel workbook

| Sheet | Contents |
|---|---|
| README | What every column means |
| Likely damaging | Every HIGH / MODERATE mutant-specific variant, priority 1 first |
| Candidate genes | One row per gene |
| Summary | Variant counts per category |

Useful columns:

| Column | Meaning |
|---|---|
| `priority` | 1 (strongest candidate) to 5 (flagged); 0 = not protein-changing |
| `warnings` | Why a variant was moved to priority 5, in plain words |
| `novel` / `known_id` | `TRUE` if the variant is not in the Ensembl/dbSNP list of known zebrafish variants; otherwise its rs number |
| `allele_check` | Amino-acid check, see "Checks on each candidate" below |
| `refseq_residues` | The amino acid at this position in published zebrafish proteins (NCBI RefSeq accession:residue) |
| `ortholog_residues` | The amino acid at this position in other vertebrates, e.g. `human:A, mouse:A, medaka:P` |
| `af_difference` | Share of mutant reads with the variant minus share of control reads (1 = only in the mutant) |
| `chrom`, `pos` | Where the variant is (chromosome, position) |
| `ref` → `alt` | The reference letter(s) and the letter(s) found in your fish |
| `gene` | The gene it falls in |
| `effect` | e.g. `stop_gained`, `missense_variant`, `frameshift_variant`, `splice_donor_variant` |
| `hgvs_p` | The protein change, e.g. `p.Arg123*` (arginine 123 becomes a stop codon) |
| `<sample>_genotype` | heterozygous / homozygous / ref |
| `<sample>_alt_reads` | Number of reads supporting the variant |
| `<sample>_alt_fraction` | Fraction of reads supporting it: ~0.5 means heterozygous, ~1.0 homozygous |

{: .tip }
In Excel use **Data → Filter** to narrow the list, for example to
`effect` = `stop_gained` and `mutant_genotype` = `homozygous`.

## Checks on each candidate

Two kinds of variants often look convincing but are almost never the mutation. The
pipeline checks every protein-changing candidate for both and moves the suspects to
**priority 5**. They stay in the tables, with the reason in the `warnings` column.

**1. Known variants.** Each variant is looked up in Ensembl's list of known zebrafish
variation (dbSNP/EVA), matching both position and changed letter. Known variants are
differences between ordinary laboratory strains. A new mutation, for example from an
ENU screen, is essentially always **novel**.

**2. Amino-acid check** (missense changes only). The zebrafish protein is lined up
with two kinds of sequence, and the amino acid at the variant position is read off
each one:

* **Published zebrafish proteins (NCBI RefSeq).** These are built largely from
  cDNAs of other fish. If one of them already has the "mutant" amino acid, then the
  reference genome (GRCz11) carries a rare allele or an error at that spot, and your
  mutant is simply normal. This is common. In one screen, three missense
  "mutations" in *aldh1a2* turned out to be exactly this.
* **The same protein in six other vertebrates** (human, rat, chicken, frog, medaka
  and fugu, from Ensembl Compara).

| `allele_check` | Meaning | Priority |
|---|---|---|
| `matches_zebrafish_refseq` | Mutant amino acid is in a published zebrafish protein: reference error or common allele | moved to 5 |
| `matches_vertebrate_consensus` | Mutant amino acid is the one most vertebrates have: probably harmless | moved to 5 |
| `alt_seen_in_orthologs` | Mutant amino acid occurs in at least one other vertebrate: may be tolerated | unchanged |
| `conserved` | Reference amino acid is the same in every species compared: change more likely damaging | unchanged, listed first |
| `variable` | The position differs between vertebrates | unchanged |
| `no_orthologs` | No other species has this gene (often zebrafish-specific genes); the mutant amino acid is not in any zebrafish RefSeq protein | unchanged |
| `not_checked` | Not a single amino-acid change, or nothing could be aligned | unchanged |

{: .note }
These checks are evidence, not proof. A priority-5 variant in a gene that fits your
phenotype perfectly is still worth a look in IGV.

## How much should I trust a variant?

A variant is **more believable** when:

* it is **high confidence**: the control was covered by enough reads to show the
  variant is absent there;
* it is **novel** and its amino-acid change is `conserved` (see above);
* it has **several supporting reads** in the mutant (5 or more is reassuring);
* its alt fraction fits the genotype: ~1.0 for homozygous, ~0.5 for heterozygous;
* the gene's known biology fits your phenotype (check ZFIN);
* it isn't sitting among dozens of other variants, a sign of a hard-to-map or
  repetitive region.

**Look at it yourself in IGV.** Open [IGV](https://igv.org/app/), choose
*Genome → Zebrafish (GRCz11)*, and load `results/alignments/<sample>.bam` together
with its `.bai` file. Then type the location, e.g. `chr5:1234567`. Real variants show
the same changed letter in many reads, from both strands.

**Then confirm it in the lab**: PCR and Sanger-sequence the variant in more mutants
and siblings. A causal recessive mutation should be homozygous in every mutant and
never homozygous in wild-type siblings.

## All the files

```
results/
├── report.html                       ← start here
├── mutant_specific/
│   ├── mutant_specific_variants.xlsx ← start here (Excel)
│   ├── peak_candidates.tsv           protein-changing variants inside the homozygosity peak
│   ├── candidate_genes.tsv
│   ├── mutant_specific_variants.tsv.gz   every mutant-specific variant, incl. non-coding (large)
│   └── mutant_specific_variants.vcf.gz   same, as VCF (for IGV)
├── homozygosity/
│   └── roh_runs.tsv                      runs of homozygosity per sample (bcftools roh)
├── structural_variants/
│   └── mutant_specific_sv.tsv
├── plots/                            figures (PNG + PDF) and the per-window tables behind them
│   ├── manhattan / variant_density / homozygosity / homozygosity_runs / coverage / candidate_map / summary
│   ├── windows_1Mb.tsv                   variant counts, homozygosity and runs of homozygosity per 1 Mb
│   └── coverage_500kb_ratio.tsv          depth and mutant/control ratio per 500 kb
├── variants/
│   ├── all_samples.annotated.vcf.gz      all variants, all samples
│   ├── all_variants_classified.tsv.gz
│   └── snpeff_summary.html
├── alignments/                           BAM files (open in IGV)
├── qc/
│   ├── multiqc_report.html               detailed quality control
│   └── fastp/                            per-sample read-cleaning reports
└── logs/                                 for troubleshooting
```
