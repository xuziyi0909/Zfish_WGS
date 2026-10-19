# Read trimming/QC, alignment, duplicate marking and alignment QC.

# Illumina adapter prefix shared by read 1 and read 2 (TruSeq / Nextera-free kits).
ADAPTER = "AGATCGGAAGAGC"


rule fastp:
    input:
        r1=lambda w: SAMPLES_DF.loc[w.sample, "fastq_1"],
        r2=lambda w: SAMPLES_DF.loc[w.sample, "fastq_2"] or [],
    output:
        r1=temp(f"{OUT}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{OUT}/trimmed/{{sample}}_R2.fastq.gz"),
        html=f"{OUT}/qc/fastp/{{sample}}.fastp.html",
        json=f"{OUT}/qc/fastp/{{sample}}.fastp.json",
    params:
        min_len=config["trimming"]["min_length"],
        adapter=ADAPTER,
    threads: min(THREADS, 16)
    log:
        f"{OUT}/logs/fastp/{{sample}}.log",
    run:
        common = (
            f"--length_required {params.min_len} --trim_poly_g --thread {threads} "
            f"--html {output.html} --json {output.json} --report_title {wildcards.sample}"
        )
        if input.r2:
            shell(
                "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
                "--detect_adapter_for_pe " + common + " > {log} 2>&1"
            )
        else:
            shell(
                "fastp -i {input.r1} -o {output.r1} --adapter_sequence {params.adapter} "
                + common + " > {log} 2>&1 && : | gzip > {output.r2}"
            )


rule align:
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2.fastq.gz",
        idx=rules.bwa_index.output,
        fai=f"{GENOME}.fai",
    output:
        bam=f"{OUT}/alignments/{{sample}}.bam",
        bai=f"{OUT}/alignments/{{sample}}.bam.bai",
        dup=f"{OUT}/qc/markdup/{{sample}}.markdup.txt",
    params:
        reads=lambda w, input: f"{input.r1} {input.r2}" if is_paired(w.sample) else input.r1,
        rg=lambda w: f"@RG\\tID:{w.sample}\\tSM:{w.sample}\\tLB:{w.sample}\\tPL:ILLUMINA",
        tmp=lambda w: f"{OUT}/alignments/tmp_{w.sample}",
    threads: THREADS
    resources:
        mem_mb=40000,
    log:
        f"{OUT}/logs/align/{{sample}}.log",
    shell:
        """
        mkdir -p {params.tmp}
        (bwa-mem2 mem -t {threads} -R '{params.rg}' {GENOME} {params.reads} \
          | samtools fixmate -m -u - - \
          | samtools sort -u -@ 4 -m 2G -T {params.tmp}/sort - \
          | samtools markdup -@ 4 -f {output.dup} -T {params.tmp}/dup --write-index - {output.bam}##idx##{output.bai}
        ) > {log} 2>&1
        rm -rf {params.tmp}
        """


rule samtools_stats:
    input:
        bam=f"{OUT}/alignments/{{sample}}.bam",
        fa=GENOME,
    output:
        stats=f"{OUT}/qc/samtools/{{sample}}.stats.txt",
        flagstat=f"{OUT}/qc/samtools/{{sample}}.flagstat.txt",
        coverage=f"{OUT}/qc/samtools/{{sample}}.coverage.txt",
    threads: 4
    params:
        mapq=config["alignment"]["min_mapping_quality"],
    shell:
        """
        samtools stats -@ {threads} --reference {input.fa} {input.bam} > {output.stats}
        samtools flagstat -@ {threads} {input.bam} > {output.flagstat}
        samtools coverage -q {params.mapq} --ff UNMAP,SECONDARY,QCFAIL,DUP {input.bam} > {output.coverage}
        """


rule multiqc:
    input:
        expand(f"{OUT}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES),
        expand(f"{OUT}/qc/samtools/{{sample}}.stats.txt", sample=SAMPLES),
        expand(f"{OUT}/qc/samtools/{{sample}}.flagstat.txt", sample=SAMPLES),
        f"{OUT}/variants/snpeff_summary.csv",
    output:
        f"{OUT}/qc/multiqc_report.html",
    params:
        outdir=f"{OUT}/qc",
    log:
        f"{OUT}/logs/multiqc.log",
    shell:
        "multiqc --force --quiet -o {params.outdir} -n multiqc_report.html {input} > {log} 2>&1"
