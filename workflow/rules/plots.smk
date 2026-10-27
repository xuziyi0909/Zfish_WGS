# Genome-wide figures: Manhattan, variant density, homozygosity, runs of homozygosity, coverage,
# candidate map.

PLOTS = ["manhattan", "variant_density", "homozygosity", "coverage", "candidate_map", "summary"]
if ROH_ENABLED:
    PLOTS.insert(3, "homozygosity_runs")


rule coverage_windows:
    """Mean read depth in 500 kb windows along the main chromosomes."""
    input:
        bams=BAMS,
        bais=expand(f"{OUT}/alignments/{{sample}}.bam.bai", sample=SAMPLES),
        fai=f"{GENOME}.fai",
    output:
        bed=temp(f"{OUT}/plots/windows_500kb.bed"),
        tsv=f"{OUT}/plots/coverage_500kb.tsv",
    params:
        mapq=config["alignment"]["min_mapping_quality"],
        header="\\t".join(["chrom", "start", "end"] + SAMPLES),
    threads: 2
    shell:
        """
        awk -v OFS='\\t' '$1 != "MT" && $2 >= 10000000 {{
            for (s = 0; s < $2; s += 500000) print $1, s, (s + 500000 < $2 ? s + 500000 : $2) }}' \
            {input.fai} > {output.bed}
        (printf '{params.header}\\n'
         samtools bedcov -Q {params.mapq} -j {output.bed} {input.bams}) > {output.tsv}
        """


rule plots:
    input:
        variants=f"{OUT}/variants/all_variants_classified.tsv.gz",
        mutant=f"{OUT}/mutant_specific/mutant_specific_variants.tsv.gz",
        coverage=f"{OUT}/plots/coverage_500kb.tsv",
        fai=f"{GENOME}.fai",
        sv=f"{OUT}/structural_variants/mutant_specific_sv.tsv" if SV_ENABLED else [],
        roh=ROH_TSV if ROH_ENABLED else [],
    output:
        png=expand(f"{OUT}/plots/{{name}}.png", name=PLOTS),
        pdf=expand(f"{OUT}/plots/{{name}}.pdf", name=PLOTS),
        windows=f"{OUT}/plots/windows_1Mb.tsv",
        cov_ratio=f"{OUT}/plots/coverage_500kb_ratio.tsv",
        peak=f"{OUT}/plots/peak.json",
        peak_candidates=f"{OUT}/mutant_specific/peak_candidates.tsv",
    params:
        controls=CONTROLS,
        mutants=MUTANTS,
        outdir=f"{OUT}/plots",
    log:
        f"{OUT}/logs/plots.log",
    script:
        "../scripts/make_plots.py"
