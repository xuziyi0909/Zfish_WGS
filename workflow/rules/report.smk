# Final plain-language HTML report.

rule report:
    input:
        fastp=expand(f"{OUT}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
        flagstat=expand(f"{OUT}/qc/samtools/{{sample}}.flagstat.txt", sample=SAMPLES),
        coverage=expand(f"{OUT}/qc/samtools/{{sample}}.coverage.txt", sample=SAMPLES),
        markdup=expand(f"{OUT}/qc/markdup/{{sample}}.markdup.txt", sample=SAMPLES),
        summary=f"{OUT}/mutant_specific/summary.json",
        genes=f"{OUT}/mutant_specific/candidate_genes.tsv",
        variants=f"{OUT}/mutant_specific/mutant_specific_variants.tsv.gz",
        sv=f"{OUT}/structural_variants/mutant_specific_sv.tsv" if SV_ENABLED else [],
        plots=expand(f"{OUT}/plots/{{name}}.png", name=PLOTS),
        peak=f"{OUT}/plots/peak.json",
        peak_candidates=f"{OUT}/mutant_specific/peak_candidates.tsv",
    output:
        f"{OUT}/report.html",
    params:
        samples=SAMPLES_DF.to_dict("records"),
        config=config,
    log:
        f"{OUT}/logs/report.log",
    script:
        "../scripts/make_report.py"
