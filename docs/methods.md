---
title: Methods & citations
nav_order: 9
---

# Methods text for your paper

You can adapt this paragraph for your methods section. The numbers in brackets are
the default settings; if you changed any, use your values from `config/config.yaml`.

> Whole-genome sequencing reads were processed with a Snakemake (v7.32) workflow.
> Adapters, poly-G tails and low-quality bases were removed with fastp (v0.23.4),
> discarding reads shorter than [50] bp. Reads were aligned to the zebrafish reference
> genome GRCz11 (Ensembl release 113) with bwa-mem2 (v2.0), and PCR duplicates were
> marked with samtools markdup (v1.19). Single-nucleotide variants and small indels
> were called jointly across all samples with bcftools mpileup/call (v1.19; mapping
> quality ≥ [20], base quality ≥ [20]), left-normalised, split into biallelic records,
> and filtered to QUAL ≥ [30]. Variant effects were predicted with snpEff (v5.2) using
> the Ensembl 113 gene annotation. A variant was considered mutant-specific if the
> mutant carried a non-reference genotype supported by ≥ [2] reads and every control
> had no supporting reads at a depth of ≥ [4] reads.
> *[Sibling mode instead:]* Candidates had ≥ [2] supporting reads and an
> alternate-allele fraction ≥ [0.9] in the mutant, and an alternate-allele fraction
> ≤ [0.6] in each sibling control covered by ≥ [4] reads.
> Variants were matched by position and allele to known zebrafish variation from
> Ensembl 113 (dbSNP/EVA). For missense variants, the zebrafish protein was aligned
> (Biopython PairwiseAligner, BLOSUM62) to the NCBI RefSeq proteins of the same gene
> and to Ensembl Compara orthologs from human, rat, chicken, *Xenopus tropicalis*,
> medaka and fugu. Variants that were known, or whose mutant residue was present in
> a zebrafish RefSeq protein or in the majority of orthologs, were deprioritised.
> Structural variants were called
> with DELLY and retained when they passed filters, were non-reference in the mutant
> and homozygous-reference in all controls.

## Please cite the tools

* **fastp**: Chen S, et al. *Bioinformatics* 2018;34:i884–i890.
* **bwa-mem2**: Vasimuddin M, et al. *IEEE IPDPS* 2019.
* **SAMtools / BCFtools**: Danecek P, et al. *GigaScience* 2021;10:giab008.
* **snpEff**: Cingolani P, et al. *Fly* 2012;6:80–92.
* **DELLY**: Rausch T, et al. *Bioinformatics* 2012;28:i333–i339.
* **MultiQC**: Ewels P, et al. *Bioinformatics* 2016;32:3047–3048.
* **Biopython**: Cock PJA, et al. *Bioinformatics* 2009;25:1422–1423.
* **Ensembl Variation / Compara**: Harrison PW, et al. *Nucleic Acids Res* 2024;52:D891–D899 (Ensembl 2024).
* **NCBI RefSeq**: O'Leary NA, et al. *Nucleic Acids Res* 2016;44:D733–D745.
* **Snakemake**: Mölder F, et al. *F1000Research* 2021;10:33.
* **Reference genome**: Howe K, et al. *Nature* 2013;496:498–503 (GRCz11, Ensembl).
