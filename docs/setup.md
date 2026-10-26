---
title: First-time setup
nav_order: 3
---

# First-time setup
{: .no_toc }

You do this **once**. It takes about 5 minutes.

1. TOC
{:toc}

---

## 1. Log in to the lab server

**Mac:** open the **Terminal** app (press ⌘+Space and type "Terminal").
**Windows:** open **PowerShell** from the Start menu.

Then type the following, using your own username, and press Enter:

```bash
ssh YOUR_USERNAME@SERVER_ADDRESS
```

Type your password when asked. The screen won't show it as you type; that's normal.

{: .note }
If you're not sure of your username or the server address, ask the lab manager.

## 2. Get your own copy of the pipeline

```bash
cd /stor/work/Gray/projects
git clone https://github.com/xuziyi0909/Zfish_WGS.git my_wgs_project
cd my_wgs_project
```

Replace `my_wgs_project` with any name you like, without spaces, for example
`jaws_mutant_2026`. Each experiment should get its own copy.

## 3. Check that everything is installed

```bash
./run_wgs.sh --check
```

Every line should say `ok`. If you see `MISSING`, send the output to the server
administrator. Here is the software the pipeline uses:

| Program | Used for |
|---|---|
| fastp ≥ 0.23 | cleaning reads |
| bwa-mem2 ≥ 2.0 | aligning reads |
| samtools / bcftools ≥ 1.17 | handling alignments, calling variants |
| snpEff ≥ 5.1 | predicting variant effects |
| DELLY ≥ 1.1 | structural variants |
| MultiQC | quality-control report |
| Snakemake 7.x, Python 3 with pandas, PyYAML, openpyxl, SciPy and Matplotlib | running the steps in order, tables and plots |
| Biopython | checking amino-acid changes against other species |

{: .tip }
**For administrators.** On Ubuntu all of these are available through `apt` (or
`conda install -c bioconda`). The pipeline creates a small private Python
environment (`.venv`) in its own folder to supply `setuptools<81`, which the
Ubuntu Snakemake package needs. Nothing is installed system-wide.

## 4. The reference genome

The first run downloads the zebrafish genome (GRCz11, Ensembl 113) and prepares it,
which takes about 30 minutes. It also downloads the list of known zebrafish variants
(about 200 MB) and protein sequences of zebrafish and six other vertebrates (about
1 GB to download, 100 MB kept), used to judge each candidate. Later runs in the **same folder** reuse it.

To avoid repeating this download for a new project, copy the prepared `reference`
folder from a previous project:

```bash
cp -r /stor/work/Gray/projects/Zfish_WGS/reference ./
find reference -type f -exec touch -d "$(date -R)" {} +
sed -i "s|^data.dir = .*|data.dir = $PWD/reference/snpeff|" reference/snpeff/snpEff.config
```

The copy takes a few minutes. Run all three lines. The second line marks the copied
files as up to date, otherwise the pipeline sees them as older than its own scripts
and rebuilds the genome index (about an hour). The third points the annotation tool at
your copy instead of the original folder.
