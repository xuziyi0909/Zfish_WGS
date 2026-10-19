---
title: Step-by-step guide
nav_order: 4
---

# Step-by-step: analyse your own data
{: .no_toc }

1. TOC
{:toc}

---

## Step 1 — Put your sequencing files on the server

Copy your `.fastq.gz` files into the `data/raw` folder inside your pipeline folder.

**From your own computer (Mac or Windows):** use a file-transfer program such as
[Cyberduck](https://cyberduck.io/) or [FileZilla](https://filezilla-project.org/).
Connect with *SFTP* to the server using your normal login, then drag the files into
`.../my_wgs_project/data/raw/`.

**If the files are already on the server** (for example in a zip from the core):

```bash
mkdir -p data/raw
unzip "/path/to/your file.zip" -d data/raw
```

{: .warning }
Don't rename the files in a way that loses `R1` and `R2`. You need to know which
file is read 1 and which is read 2.

## Step 2 — Tell the pipeline about your samples

Open the sample sheet:

```bash
nano config/samples.tsv
```

It has **four columns separated by TAB characters** (not spaces):

| sample | group | fastq_1 | fastq_2 |
|---|---|---|---|
| a short name, no spaces | `control` or `mutant` | path to the read-1 file | path to the read-2 file (leave empty for single-end) |

**Example: paired-end, two siblings and two mutants**

```
sample	group	fastq_1	fastq_2
wt_sib1	control	data/raw/WT1_S1_L001_R1_001.fastq.gz	data/raw/WT1_S1_L001_R2_001.fastq.gz
wt_sib2	control	data/raw/WT2_S2_L001_R1_001.fastq.gz	data/raw/WT2_S2_L001_R2_001.fastq.gz
mut1	mutant	data/raw/MUT1_S3_L001_R1_001.fastq.gz	data/raw/MUT1_S3_L001_R2_001.fastq.gz
mut2	mutant	data/raw/MUT2_S4_L001_R1_001.fastq.gz	data/raw/MUT2_S4_L001_R2_001.fastq.gz
```

**Example: single-end (one file per sample)**

```
sample	group	fastq_1	fastq_2
control_0038	control	data/raw/L5756-0013-0038_S33_L005_R2_001.fastq.gz	
mutant_0039	mutant	data/raw/L5756-0013-0039_S34_L005_R1_001.fastq.gz	
```

In **nano**: move with the arrow keys, press the **Tab** key between columns, save with
**Ctrl+O** then **Enter**, and exit with **Ctrl+X**.

{: .tip }
**Prefer Excel?** Fill in the four columns in Excel, choose *File → Save As →
"Tab delimited text (.txt)"*, then upload the file as `config/samples.tsv`.

{: .note }
**One sample sequenced on several lanes** (files with `L001`, `L002`, ...)? Join them
first, keeping R1 and R2 separate:
`cat S1_L001_R1_001.fastq.gz S1_L002_R1_001.fastq.gz > S1_R1.fastq.gz`

Now check the sheet:

```bash
./run_wgs.sh --check
```

You should see one `ok` line per sample. Otherwise the message tells you what to fix,
for example `file not found` or `group must be 'control' or 'mutant'`.

## Step 3 — Run it

The analysis takes hours, so start it inside **screen**. That way it keeps running
after you close your laptop.

```bash
screen -S wgs          # open a "screen" session named wgs
./run_wgs.sh           # start the pipeline
```

Now press **Ctrl+A**, then **D**, to leave it running in the background. You can
log out.

To check on it later:

```bash
screen -r wgs          # re-attach and see the progress
```

When it finishes you will see:

```
Finished! Open this file in a web browser:  results/report.html
```

{: .tip }
If it stopped, for example after a power cut or with an error you have since fixed,
run `./run_wgs.sh` again. Finished steps are skipped and it continues where it
stopped.

### How long does it take?

Measured on the lab server (72 cores, `--cores 48`) for the example dataset of
two single-end samples, 61 M and 74 M reads:

| Step | Time |
|---|---|
| Download and index the genome (first time only) | ~25 min |
| Cleaning and aligning (both samples in parallel) | ~45 min |
| Variant calling | ~10 min |
| Structural variants (runs alongside variant calling) | ~25 min |
| Annotation, comparison and reports | ~35 min |
| **Total** | **about 2 hours** (about 1.5 h when the genome is already prepared) |

Deeper or paired-end sequencing, or more samples, takes proportionally longer:
plan on about 1 hour of alignment per 100 M reads.

## Step 4 — Look at your results

Download the `results` folder to your computer using Cyberduck or FileZilla, then
double-click:

1. **`results/report.html`**: the summary report. It opens in any web browser.
2. **`results/mutant_specific/mutant_specific_variants.xlsx`**: your candidate
   mutations, in Excel.

[How to read the results →](results.html)

## Optional: change the settings

The defaults suit most experiments. To change them, open `config/config.yaml`; every
setting is explained there. See [Settings](settings.html).
