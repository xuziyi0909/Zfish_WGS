# Joint small-variant calling, annotation, mutant-vs-control comparison and SVs.

BAMS = expand(f"{OUT}/alignments/{{sample}}.bam", sample=SAMPLES)
ANNOTATED_VCF = f"{OUT}/variants/all_samples.annotated.vcf.gz"
# with known-variant annotation, snpEff writes an intermediate file first
SNPEFF_VCF = f"{OUT}/variants/all_samples.snpeff.vcf.gz" if KNOWN_ENABLED else ANNOTATED_VCF


rule call_chunk:
    input:
        bams=BAMS,
        bais=expand(f"{OUT}/alignments/{{sample}}.bam.bai", sample=SAMPLES),
        bed=f"{REF_DIR}/chunks/chunk_{{chunk}}.bed",
        fai=f"{GENOME}.fai",
    output:
        temp(f"{OUT}/variants/chunks/chunk_{{chunk}}.bcf"),
    params:
        mapq=config["alignment"]["min_mapping_quality"],
        baseq=config["variant_calling"]["min_base_quality"],
    log:
        f"{OUT}/logs/call/chunk_{{chunk}}.log",
    shell:
        """
        (bcftools mpileup -Ou -f {GENOME} -R {input.bed} -q {params.mapq} -Q {params.baseq} \
            --max-depth 1000 -a FORMAT/AD,FORMAT/DP {input.bams} \
          | bcftools call -mv -Ob -o {output}) > {log} 2>&1
        """


rule merge_calls:
    input:
        expand(f"{OUT}/variants/chunks/chunk_{{chunk}}.bcf", chunk=CHUNKS),
    output:
        vcf=f"{OUT}/variants/all_samples.raw.vcf.gz",
    threads: 4
    log:
        f"{OUT}/logs/call/merge.log",
    shell:
        """
        (bcftools concat --threads {threads} -Ou {input} \
          | bcftools norm --threads {threads} -f {GENOME} -m -any -Ou \
          | bcftools norm --threads {threads} -d exact -Oz -o {output.vcf} \
          && bcftools index -t {output.vcf}) > {log} 2>&1
        """


rule filter_calls:
    input:
        f"{OUT}/variants/all_samples.raw.vcf.gz",
    output:
        f"{OUT}/variants/all_samples.filtered.vcf.gz",
    params:
        qual=config["variant_calling"]["min_qual"],
    shell:
        """
        bcftools view -i 'QUAL>={params.qual} && ALT!="*"' -Oz -o {output} {input}
        bcftools index -t {output}
        """


rule homozygosity_runs:
    """Runs of homozygosity per sample (bcftools roh), on SNP genotypes with weak calls removed."""
    input:
        f"{OUT}/variants/all_samples.filtered.vcf.gz",
    output:
        ROH_TSV,
    params:
        min_reads=config.get("homozygosity_runs", {}).get("min_reads", 3),
        min_allele=config.get("homozygosity_runs", {}).get("min_reads_each_allele", 2),
        af=config.get("homozygosity_runs", {}).get("allele_frequency", 0.4),
        gt_error=config.get("homozygosity_runs", {}).get("genotype_error_phred", 30),
        az_to_hw=config.get("homozygosity_runs", {}).get("az_to_hw", 1e-10),
    log:
        f"{OUT}/logs/homozygosity_runs.log",
    shell:
        # setGT: a call with too few reads, or a "heterozygote" resting on a single read of one
        # allele (usually an error-prone or misplaced read), is set to missing before the HMM
        """
        (bcftools view -v snps -Ou {input} \
          | bcftools +setGT -Ou - -- -t q -n . \
              -i 'FMT/DP<{params.min_reads} | (GT="het" & (FMT/AD[:0]<{params.min_allele} | FMT/AD[:1]<{params.min_allele}))' \
          | bcftools roh --AF-dflt {params.af} -G {params.gt_error} --ignore-homref -H {params.az_to_hw} -O r - \
          | awk -v OFS='\\t' 'BEGIN {{print "sample", "chrom", "start", "end", "length", "sites", "quality"}}
                              $1 == "RG" {{print $2, $3, $4, $5, $6, $7, $8}}' > {output}) 2> {log}
        """


rule annotate:
    input:
        vcf=f"{OUT}/variants/all_samples.filtered.vcf.gz",
        db=SNPEFF_DB,
        config=SNPEFF_CONFIG,
    output:
        vcf=temp(SNPEFF_VCF) if KNOWN_ENABLED else SNPEFF_VCF,
        tbi=temp(f"{SNPEFF_VCF}.tbi") if KNOWN_ENABLED else f"{SNPEFF_VCF}.tbi",
        csv=f"{OUT}/variants/snpeff_summary.csv",
        html=f"{OUT}/variants/snpeff_summary.html",
    params:
        name=REF_NAME,
    log:
        f"{OUT}/logs/snpeff.log",
    shell:
        """
        (snpEff -Xmx16g ann -c {input.config} -nodownload -csvStats {output.csv} \
            -stats {output.html} {params.name} {input.vcf} \
          | sed -n '/^##fileformat/,$p' \
          | bgzip > {output.vcf} && tabix -p vcf {output.vcf}) 2> {log}
        """


rule annotate_known:
    """Fill the ID column with the rsID of known Ensembl/dbSNP variants (same position and allele)."""
    input:
        vcf=SNPEFF_VCF,
        tbi=f"{SNPEFF_VCF}.tbi",  # bcftools annotate needs an indexed target
        known=KNOWN_VCF,
        known_tbi=f"{KNOWN_VCF}.tbi",
    output:
        vcf=ANNOTATED_VCF,
        tbi=f"{ANNOTATED_VCF}.tbi",
    threads: 4
    log:
        f"{OUT}/logs/annotate_known.log",
    shell:
        """
        (bcftools annotate --threads {threads} -a {input.known} -c ID --pair-logic exact \
            -Oz -o {output.vcf} {input.vcf} && tabix -f -p vcf {output.vcf}) > {log} 2>&1
        """


rule compare_samples:
    """Classify every variant as mutant-specific / shared / control-only."""
    input:
        vcf=ANNOTATED_VCF,
        protein_refs=list(PROTEIN_REFS.values()) if PROTEIN_CHECK else [],
    output:
        all_tsv=f"{OUT}/variants/all_variants_classified.tsv.gz",
        mut_tsv=f"{OUT}/mutant_specific/mutant_specific_variants.tsv.gz",
        genes=f"{OUT}/mutant_specific/candidate_genes.tsv",
        xlsx=f"{OUT}/mutant_specific/mutant_specific_variants.xlsx",
        vcf=f"{OUT}/mutant_specific/mutant_specific_variants.vcf.gz",
        summary=f"{OUT}/mutant_specific/summary.json",
    params:
        controls=CONTROLS,
        mutants=MUTANTS,
        min_alt=config["comparison"]["min_alt_reads_mutant"],
        min_ctrl_depth=config["comparison"]["min_depth_control"],
        min_depth=config["variant_calling"]["min_depth"],
        rule=config["comparison"]["mutant_rule"],
        mode=COMPARISON_MODE,
        min_af_mutant=config["comparison"].get("min_alt_fraction_mutant", 0.9),
        max_af_control=config["comparison"].get("max_alt_fraction_control", 0.6),
        demote_flagged=config["comparison"].get("demote_flagged", True),
        known=KNOWN_ENABLED,
        protein_refs=PROTEIN_REFS if PROTEIN_CHECK else {},
        species=config.get("protein_check", {}).get("species", []),
    threads: 8
    log:
        f"{OUT}/logs/compare_samples.log",
    script:
        "../scripts/compare_samples.py"


rule structural_variants:
    input:
        bams=BAMS,
        fai=f"{GENOME}.fai",
    output:
        bcf=f"{OUT}/structural_variants/all_samples.sv.bcf",
    threads: 2
    log:
        f"{OUT}/logs/delly_call.log",
    shell:
        "delly call -g {GENOME} -q 20 -o {output.bcf} {input.bams} > {log} 2>&1"


rule structural_variants_mutant:
    input:
        bcf=f"{OUT}/structural_variants/all_samples.sv.bcf",
    output:
        tsv=f"{OUT}/structural_variants/mutant_specific_sv.tsv",
        vcf=f"{OUT}/structural_variants/mutant_specific_sv.vcf.gz",
    params:
        controls=CONTROLS,
        mutants=MUTANTS,
    log:
        f"{OUT}/logs/delly_filter.log",
    script:
        "../scripts/sv_compare.py"
