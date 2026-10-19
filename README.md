# Zebrafish WGS Pipeline

Find DNA changes (SNVs, indels, structural variants) present in **mutant** zebrafish
but absent from **control** fish, starting from raw `.fastq.gz` files. The pipeline
is built for people with no coding or WGS experience.

📖 **Full step-by-step guide: https://xuziyi0909.github.io/Zfish_WGS/**

## Quick start

```bash
git clone https://github.com/xuziyi0909/Zfish_WGS.git my_project && cd my_project
# 1. copy your FASTQ files into data/raw/
# 2. list them in config/samples.tsv (sample, group = control|mutant, fastq_1, fastq_2)
./run_wgs.sh --check     # checks software and the sample sheet
./run_wgs.sh             # runs everything (resumes automatically if interrupted)
```

Outputs:

- `results/report.html`: a plain-language summary (QC, candidate genes, top variants)
- `results/mutant_specific/mutant_specific_variants.xlsx`: candidate mutations for Excel

## What it does

| Step | Tool |
|---|---|
| Adapter / quality trimming | fastp |
| Alignment to GRCz11 (Ensembl 113, downloaded automatically) | bwa-mem2 |
| Duplicate marking, alignment QC | samtools |
| Joint SNV/indel calling | bcftools |
| Effect prediction | snpEff (database built locally from the Ensembl GTF) |
| Known-variant annotation (Ensembl/dbSNP rsIDs) | bcftools annotate |
| Mutant-vs-control classification (`absent` or `sibling` mode) | `workflow/scripts/compare_samples.py` |
| Missense check vs zebrafish RefSeq proteins and vertebrate orthologs | `workflow/scripts/protein_check.py` (Biopython) |
| Structural variants | DELLY |
| Reports | MultiQC + custom HTML report |

It supports single-end and paired-end data, any number of control and mutant
samples, and resumes automatically after an interruption.

## Layout

```
run_wgs.sh                 one-command launcher (checks, then runs Snakemake)
config/config.yaml         settings (commented)
config/samples.tsv         sample sheet
workflow/Snakefile         pipeline definition
workflow/rules/*.smk       pipeline steps
workflow/scripts/          helper scripts (comparison, report, SVs)
docs/                      GitHub Pages user guide
```
