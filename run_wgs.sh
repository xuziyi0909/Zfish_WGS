#!/usr/bin/env bash
# ======================================================================
#  Zebrafish WGS pipeline - one command to run everything.
#
#  Usage:
#     ./run_wgs.sh              check your setup, then run the pipeline
#     ./run_wgs.sh --check      only check your setup and sample sheet
#     ./run_wgs.sh --cores 32   use 32 CPU cores (default: half the machine)
#
#  Safe to re-run: finished steps are skipped and it picks up where it stopped.
# ======================================================================
set -uo pipefail

cd "$(dirname "$(readlink -f "$0")")"

CORES=$(( $(nproc) / 2 ))
CHECK_ONLY=0
EXTRA=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --cores) CORES="$2"; shift 2 ;;
        --check) CHECK_ONLY=1; shift ;;
        -h|--help) sed -n 2,11p "$0"; exit 0 ;;
        *) EXTRA+=("$1"); shift ;;
    esac
done

green() { printf '\033[32m%s\033[0m\n' "$*"; }
red()   { printf '\033[31m%s\033[0m\n' "$*"; }
bold()  { printf '\033[1m%s\033[0m\n' "$*"; }

bold "== Step 1/3: checking software"
missing=0
for tool in snakemake fastp bwa-mem2 samtools bcftools bgzip tabix snpEff delly multiqc python3 java wget; do
    if command -v "$tool" >/dev/null 2>&1; then echo "  ok   $tool"; else red "  MISSING  $tool"; missing=1; fi
done
for module in pandas yaml openpyxl scipy matplotlib Bio; do
    if python3 -c "import $module" >/dev/null 2>&1; then echo "  ok   python module $module"; else red "  MISSING  python module $module"; missing=1; fi
done
if [[ $missing -eq 1 ]]; then
    red "Some programs are missing. Ask your server administrator to install them (see docs/setup.md)."
    exit 1
fi

bold "== Step 2/3: checking the sample sheet"
python3 - <<'EOF' || exit 1
import os, sys, yaml, pandas as pd
cfg = yaml.safe_load(open("config/config.yaml"))
path = cfg["samples"]
if not os.path.exists(path):
    sys.exit(f"  ERROR: sample sheet '{path}' not found")
df = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")
df.columns = [c.strip() for c in df.columns]
problems = []
for col in ("sample", "group", "fastq_1"):
    if col not in df.columns:
        problems.append(f"missing column '{col}' (columns must be separated by TABs, not spaces)")
if not problems:
    if df["sample"].duplicated().any():
        problems.append("sample names must be unique: " + ", ".join(df.loc[df['sample'].duplicated(), 'sample']))
    for _, r in df.iterrows():
        name = r["sample"].strip()
        if not name or any(ch in name for ch in " /\\|"):
            problems.append(f"sample name '{name}' must not be empty or contain spaces or slashes")
        if r["group"].strip().lower() not in ("control", "mutant"):
            problems.append(f"{name}: group must be 'control' or 'mutant' (found '{r['group']}')")
        for col in ("fastq_1", "fastq_2"):
            f = r.get(col, "").strip()
            if f and not os.path.exists(f):
                problems.append(f"{name}: file not found: {f}")
    groups = set(df["group"].str.strip().str.lower())
    if "mutant" not in groups:
        problems.append("at least one sample must have group 'mutant'")
    if "control" not in groups:
        print("  WARNING: no control sample - every variant in the mutant will be reported")
if problems:
    print("  Problems found in", path)
    for p in problems:
        print("   -", p)
    sys.exit(1)
for _, r in df.iterrows():
    mode = "paired-end" if r.get("fastq_2", "").strip() else "single-end"
    print(f"  ok   {r['sample']:<20} {r['group']:<8} {mode}")
EOF

if [[ $CHECK_ONLY -eq 1 ]]; then
    green "Setup looks good. Run ./run_wgs.sh to start the analysis."
    exit 0
fi

OUT=$(python3 -c "import yaml;print(yaml.safe_load(open('config/config.yaml'))['output_dir'])")
mkdir -p "$OUT/logs"
LOG="$OUT/logs/run_$(date +%Y%m%d_%H%M%S).log"

bold "== Step 3/3: running the analysis on $CORES cores"
echo "  The first run also downloads and indexes the zebrafish genome (~1-2 hours, once)."
echo "  A full run takes several hours. You can close this window if you started it"
echo "  inside 'screen' or 'tmux' (see docs). Progress log: $LOG"

# Snakemake on Ubuntu needs pkg_resources (old setuptools); keep it in a private
# environment next to the pipeline so nothing system-wide is changed.
if [[ ! -x .venv/bin/python ]] || ! .venv/bin/python -c "import pkg_resources" 2>/dev/null; then
    echo "  Preparing a private Python environment (one time only)..."
    python3 -m venv --system-site-packages .venv && .venv/bin/pip install -q "setuptools<81" || {
        red "Could not prepare the Python environment (see docs/troubleshooting.md)."; exit 1; }
fi
export PYTHONWARNINGS=ignore

.venv/bin/python -m snakemake --snakefile workflow/Snakefile --cores "$CORES" --resources mem_mb=200000 \
    --rerun-incomplete --keep-going --printshellcmds "${EXTRA[@]}" 2>&1 | tee "$LOG"
status=${PIPESTATUS[0]}

if [[ $status -eq 0 ]]; then
    green "Finished! Open this file in a web browser:  $OUT/report.html"
    green "Main results table (Excel):                 $OUT/mutant_specific/mutant_specific_variants.xlsx"
else
    red "The pipeline stopped with an error. Look for 'Error in rule' in $LOG"
    red "and the log file it mentions. See docs/troubleshooting.md for common fixes."
fi
exit $status
