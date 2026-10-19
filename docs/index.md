---
title: Home
layout: home
nav_order: 1
---

# Find the DNA changes in your mutant zebrafish
{: .fs-8 }

A ready-to-use whole-genome sequencing (WGS) pipeline. You give it the sequencing
files for your **control** and **mutant** fish; it gives you back an **Excel table
of candidate mutations** and a **report you open in your web browser**.
{: .fs-5 .fw-300 }

[Start the step-by-step guide](run.html){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 }
[What is WGS?](concepts.html){: .btn .fs-5 .mb-4 .mb-md-0 }

---

## Do I need to know how to code?

No. You will type **three short commands**, and every one is written out for you
to copy and paste. Everything else runs on its own.

## What you need

| | |
|---|---|
| Sequencing files | `.fastq.gz` files from the sequencing core, for at least one mutant and (ideally) one control |
| An account on the lab server | All the software is already installed there |
| Time | ~15 minutes of your time; the computer then works for about 2–6 hours |

## What the pipeline does

```
 your .fastq.gz files
        │
   1. Clean reads ............ removes adapters and low-quality bases   (fastp)
   2. Map to the genome ...... places each read on the zebrafish genome (bwa-mem2)
   3. Remove duplicates ...... ignores PCR copies of the same fragment   (samtools)
   4. Find variants .......... spots differences from the reference     (bcftools)
   5. Predict effects ........ which gene? does it change the protein?  (snpEff)
   6. Compare mutant vs control  keeps changes found only in the mutant
   7. Large rearrangements ... deletions, inversions, duplications      (DELLY)
        │
   report.html  +  mutant_specific_variants.xlsx
```

## The three commands

```bash
cd /stor/work/Gray/projects/Zfish_WGS      # 1. go to the pipeline folder
nano config/samples.tsv                    # 2. list your samples
./run_wgs.sh                               # 3. run everything
```

See the [step-by-step guide](run.html) for the details.

## Example results

The pipeline was set up and tested on a real dataset: one control fish and one
mutant fish. [See what the results look like](example.html).
