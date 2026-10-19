---
title: WGS in plain words
nav_order: 2
---

# Whole-genome sequencing in plain words
{: .no_toc }

1. TOC
{:toc}

---

## The big picture

Every cell holds the fish's **genome**: about 1.4 billion DNA letters (A, C, G, T)
across 25 chromosomes. A sequencing machine cannot read a whole chromosome in one
go. Instead it reads millions of short, random pieces of about 150 letters each,
called **reads**.

The pipeline lines those pieces up against the published zebrafish genome, the
**reference**, and looks for places where your fish differs from it.

## Words you will see

| Word | Meaning |
|---|---|
| **FASTQ** (`.fastq.gz`) | The file from the sequencer. It holds every read plus a quality score for each letter. |
| **Read 1 / Read 2** (`R1`, `R2`) | With *paired-end* sequencing each DNA fragment is read from both ends, giving two files per sample. *Single-end* data has one file. The pipeline handles both. |
| **Reference genome** | The standard zebrafish genome sequence (GRCz11). The pipeline downloads it for you. |
| **Alignment / mapping** | Finding where each read belongs on the reference. |
| **BAM** | The file of aligned reads. You can view it in [IGV](https://igv.org/app/). |
| **Coverage (depth)** | How many reads cover a given DNA letter. 10× means each letter was read about 10 times. More is better: ≥15× is good, ≥8× is usable. |
| **Variant** | A position where your fish's DNA differs from the reference. |
| **SNV** | Single-nucleotide variant: one letter changed (e.g. A→G). |
| **Indel** | A few letters inserted or deleted. |
| **Structural variant (SV)** | A large change: a deletion, duplication, inversion or translocation of hundreds to millions of letters. |
| **Heterozygous** | The variant is on one of the two chromosome copies. |
| **Homozygous** | The variant is on both copies. Recessive mutant phenotypes usually need this. |
| **VCF** (`.vcf.gz`) | The standard file format for variants. |
| **Impact** | snpEff's prediction of how much a variant disrupts a gene: **HIGH** (stop codon, frameshift, splice site), **MODERATE** (amino-acid change), **LOW** (silent change), **MODIFIER** (outside a protein-coding region). |

## Why compare against a control?

Any two zebrafish differ at **millions** of positions, even siblings from the same
lab line. Almost all of these differences are harmless background. Variants that
also appear in the control cannot explain the mutant phenotype, so the pipeline
removes them. What remains is a much shorter list of **mutant-specific** variants.

{: .tip }
Your best control is a **wild-type sibling** (or a pool of siblings) from the same
cross as the mutant. The more closely related the control, the less background
is left over. With sibling controls, set `mode: sibling` in the settings
([why](settings.html#choose-a-comparison-mode)): a wild-type-looking sibling can carry
one copy of a recessive mutation.

## What the pipeline cannot do for you

* It cannot prove that a variant **causes** the phenotype. Confirm your top
  candidates by PCR and Sanger sequencing, by checking that they segregate with
  the phenotype, or with a second allele or CRISPR knockout.
* At low coverage (below about 10×) it can miss one copy of a heterozygous variant
  in the control. Some "mutant-specific" calls are therefore really shared.
  [More on reading the results](results.html#how-much-should-i-trust-a-variant).
