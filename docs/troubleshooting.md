---
title: Troubleshooting & FAQ
nav_order: 8
---

# Troubleshooting & FAQ
{: .no_toc }

1. TOC
{:toc}

---

## `./run_wgs.sh: Permission denied`

Make the script runnable, once:

```bash
chmod +x run_wgs.sh
```

## `MISSING <program>`

A program isn't installed on the server. Send the output of `./run_wgs.sh --check`
to the server administrator; [First-time setup](setup.html) lists what's needed.

## `missing column 'group'` or `file not found`

* The columns in `config/samples.tsv` must be separated by **Tab**, not spaces.
  In nano, press the Tab key between columns.
* File paths are relative to the pipeline folder, e.g. `data/raw/name.fastq.gz`.
  Check the exact file name with `ls data/raw`.
* Spaces in file names cause trouble. Rename the file, e.g.
  `mv "data/raw/my file.fastq.gz" data/raw/my_file.fastq.gz`.

## `The pipeline stopped with an error`

1. Open the progress log it names and search for `Error in rule`.
2. Under it, a line starting with `log:` names the detailed log for the failed step.
   Open that file: the last lines usually explain the problem.
3. Fix it, then run `./run_wgs.sh` again. Finished steps are not repeated.

Common causes:

| Message in the log | Fix |
|---|---|
| `No space left on device` | Free up disk space. Each sample needs about 3× the size of its FASTQ files. |
| `gzip: unexpected end of file` / `invalid compressed data` | The FASTQ file is damaged or incomplete. Copy it again. |
| `LockException` / `Directory cannot be locked` | Another run is still going, or one crashed. If you're sure none is running: `.venv/bin/python -m snakemake --snakefile workflow/Snakefile --unlock` |
| `Killed` / `MemoryError` | Lower `threads_per_job` in `config/config.yaml`, or run with `./run_wgs.sh --cores 16` |

## I closed my laptop and the run stopped

Start it inside `screen` ([step 3](run.html#step-3--run-it)). Then simply re-run
`./run_wgs.sh`; it continues where it stopped.

## How do I start over completely?

Delete the results folder; keep `reference/` so the genome isn't downloaded again:

```bash
rm -rf results
./run_wgs.sh
```

## I have hundreds of mutant-specific variants. Which is the one?

That's normal, especially with one fish per group or low coverage. To narrow it down:

1. Start with the **Candidate genes** sheet: HIGH-impact, high-confidence,
   homozygous hits first.
2. Check each gene on ZFIN for expression or phenotypes that fit yours.
3. Look at the top hits in IGV ([how](results.html#how-much-should-i-trust-a-variant)).
4. Confirm by PCR and Sanger sequencing in more mutants and siblings.
5. For the next experiment, sequence **pools** of mutants and of siblings (20–50
   fish each) or several individuals. Background variation then drops sharply.

## My data are paired-end. Do I need to do anything?

Put the R1 file in `fastq_1` and the R2 file in `fastq_2`. The pipeline handles
the rest. Paired-end data also makes structural-variant detection far more sensitive.

## Can I use another species or genome build?

Yes. Change the `reference` section of `config/config.yaml` to point at another
Ensembl FASTA and GTF, and give it a new `name`. Then delete the `reference/`
folder, or set a new `reference.dir`.
