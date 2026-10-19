---
title: Example results
nav_order: 5
---

# Example: the dataset used to build this pipeline
{: .no_toc }

1. TOC
{:toc}

---

## The data

| Sample | Group | File | Reads |
|---|---|---|---|
| `control_0038` | control (sibling) | `L5756-0013-0038_S33_L005_R2_001.fastq.gz` | 61.3 M × 151 bp, single-end |
| `mutant_0039` | mutant | `L5756-0013-0039_S34_L005_R1_001.fastq.gz` | 73.5 M × 151 bp, single-end |

Only one read file per sample was available (R2 for the control and R1 for the
mutant), so the data was analysed as **single-end**.

The sample sheet used:

```
sample	group	fastq_1	fastq_2
control_0038	control	data/raw/L5756-0013-0038_S33_L005_R2_001.fastq.gz	
mutant_0039	mutant	data/raw/L5756-0013-0039_S34_L005_R1_001.fastq.gz	
```

## Data quality

| | control_0038 | mutant_0039 | Verdict |
|---|---|---|---|
| Reads kept after cleaning | 98.8% | 99.1% | ✓ Good |
| High-quality bases (Q30) | 97.7% | 97.3% | ✓ Good |
| Mapped to GRCz11 | 98.4% | 98.4% | ✓ Good |
| Duplicates | 13.2% | 12.8% | ✓ Good |
| **Mean usable coverage** | **2.8×** | **3.4×** | ✕ Low |
| Genome covered by ≥ 1 read | 78.6% | 81.0% | ✕ Low |

The sequencing itself is excellent; there just isn't much of it. About 89% of reads
contained adapter sequence, meaning the DNA fragments were short, often under
150 bp, so about a quarter of all sequenced letters were trimmed off. Together
with duplicates and reads that map ambiguously in repeats, that leaves about 3×
usable coverage. 15× or more is ideal.

## Variants

14.1 million variants passed the quality filter. That's normal for zebrafish: lab
strains are far more genetically diverse than mice or humans.

Fish 0038 is a wild-type-looking **sibling** of the mutant, so the run uses
`mode: sibling` ([why](settings.html#choose-a-comparison-mode)). A candidate must be
carried in at least 90% of the mutant's reads. The sibling may carry it in up to 60%
of its reads, because it may carry one copy.

| Category | Variants | HIGH impact | MODERATE impact |
|---|---:|---:|---:|
| Mutant-specific (high confidence) | 18,988 | 3 | 132 |
| Mutant-enriched (high confidence) | 154,760 | 29 | 1,028 |
| Mutant-specific (low confidence) | 2,111,050 | 379 | 10,348 |
| Mutant, not homozygous | 1,095,204 | 284 | 5,706 |
| Shared with the sibling | 5,341,182 | 1,260 | 40,493 |
| Control only | 5,104,933 | 1,175 | 22,216 |
| Low evidence | 319,791 | 104 | 765 |

The low-confidence group is large because at about 3× the sibling often has fewer
than 4 reads at a position, too few to judge.

### Checks on the candidates

Of the 1,192 high-confidence protein-changing variants:

* **623 are known** zebrafish variants (Ensembl/dbSNP): ordinary strain differences.
* **69** change the amino acid to the one already found in a published zebrafish
  protein (RefSeq). The reference genome carries a rare allele or an error there.
* **62** change it to the residue most other vertebrates have.
* **38** hit an amino acid that is identical in every species compared. These are the
  most likely to be damaging.

Known and flagged variants are moved to priority 5. Several of them looked like
excellent jaw candidates:

| Variant | Reads: mutant · sibling | Why it is priority 5 |
|---|---|---|
| *aldh1a2* P393A | 3/3 · 2/5 | Known (rs40833506); the zebrafish RefSeq protein NP_571925 has A, like the mutant |
| *slc39a13* A23P | 4/4 · 0/4 | Known polymorphism (rs507894324), despite a small jaw in the mouse knockout |
| *lgr4* I745V | 7/7 · 0/0 | V is the residue in human, rat, chicken, frog, medaka and fugu |
| *mmp14a* N338D | 3/3 · 1/1 | Known (rs179907423); the zebrafish RefSeq protein has D |
| *b3gat3* S94A | 9/9 · 3/3 | Known (rs180151941); the sibling looks homozygous too |

### Top candidates

| Priority | Variants | Examples |
|---|---:|---|
| 1 (HIGH, high confidence, homozygous) | 21 | stop gained in *SLC35F4*, *slc52a2*, *gig2f*; splice site in *itga8*, *TPR*, *cyt1l*; frameshift in *CENPV*, *CCM2L*; stop lost in *si:dkey-1h4.4* (chr7) |
| 2 (MODERATE, high confidence, homozygous) | 503 | 489 missense, 14 in-frame insertions or deletions |
| 3 | 3 | |
| 4 (low confidence) | 3,268 | |
| 5 (flagged) | 8,124 | known variants and flagged amino-acid changes |

Each priority-1 variant rests on only 3–8 reads in the mutant. Most also show the
variant in some sibling reads, as expected for a carrier, but at this depth that
can't be told apart from chance. They're **leads to test by PCR and Sanger
sequencing**, not conclusions. Only one of the 21 lies in the chr7 peak below.

A useful filter in the Excel file is priority 1–2, `novel` = TRUE and
`allele_check` = `conserved`. That leaves 22 variants genome-wide, and only two of
them inside the chr7 peak: *tmppe* R370W and *TJP1* P1240L.

### The mapping signal: chromosome 7

![Variant density](assets/variant_density.png)

The clearest result isn't any single variant but a **region**. Chromosome 7,
**22–33 Mb**, reaches up to 2,582 homozygous candidate variants per Mb, about 20
times the genome average of 125. It holds 14,867 of them. The homozygosity
difference track and the Manhattan plot peak at the same place, and the region
shows no change in copy number.

![Homozygosity](assets/homozygosity.png)

That's the expected footprint of a **recessive mutation**: the mutant carries two
identical copies of the chromosome segment around it, while the sibling doesn't.
The region holds 44 priority 1–2 variants, all novel. Most are passengers carried
along with the causal change. The leads worth testing first:

| Variant (chr7) | Position | Reads: mutant · sibling | Notes |
|---|---|---|---|
| *tmppe* R370W | 29,157,506 | 5/5 · 0/6 | Novel; R in human, rat and chicken (conserved) |
| *TJP1* P1240L | 31,126,090 | 3/3 · 2/8 | Novel; P in human, rat, chicken and frog (conserved); the sibling looks like a carrier |
| *HERC1* Q3550H | 29,211,339 | 7/7 · 0/5 | Novel, best read support; not conserved (other vertebrates have E) |
| *si:dkey-1h4.4* stop lost | 30,929,805 | 4/4 · 1/4 | The only priority-1 variant in the peak; unnamed gene |

The known jaw genes *cog4* and *thap11* sit in the core of the peak with no coding
change; *thap11* has a novel 5′ UTR variant. The full write-up is in
`analysis/jaw_candidates/STATUS.md`.

{: .note }
Genotyping a few chr7 variants (for example *HERC1* Q3550H, plus markers near 22 and
33 Mb) in about 20 mutants and 20 siblings will quickly show whether the region is
linked to the phenotype, and narrow it down.

![Manhattan plot](assets/manhattan.png)

### Structural variants

228 mutant-specific structural variants passed the filters: 214 deletions, 8
duplications, 4 insertions and 2 inversions. 148 of them overlap a gene. With
single-end reads only split reads can be used, so many larger events will be
missed.

## What this example teaches

{: .warning }
With **one fish per group at about 3× coverage**, most "mutant-specific" variants
are ordinary differences between two individuals, plus variants missed in the
thinly covered control. The pipeline flags all of this in the report.

To get a much cleaner answer from the next experiment:

1. **Sequence both reads (R1 + R2)** of each sample. That alone roughly doubles
   coverage and makes structural-variant calls far more sensitive.
2. **Aim for ≥ 15× per sample**, about 25 Gb of sequence for zebrafish.
3. **Use related controls**, as here: wild-type siblings from the same cross,
   ideally several, or a pool.
4. For a recessive mutation from a screen, **pool 20–50 mutants and 20–50 siblings**.
   The causal mutation is then the variant that is homozygous in the mutant pool
   and rarely homozygous in the sibling pool. Run such pools with `mode: sibling`
   ([settings](settings.html#choose-a-comparison-mode)).
