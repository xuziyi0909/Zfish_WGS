# Reference genome download and indexing (runs once, then reused).

rule download_reference:
    output:
        fasta=GENOME,
        gtf=GTF,
    params:
        fasta_url=config["reference"]["fasta_url"],
        gtf_url=config["reference"]["gtf_url"],
    log:
        f"{REF_DIR}/logs/download.log",
    shell:
        """
        (wget -q -c -O {output.fasta}.gz {params.fasta_url}
         wget -q -c -O {output.gtf}.gz {params.gtf_url}
         gzip -dc {output.fasta}.gz > {output.fasta}.tmp && mv {output.fasta}.tmp {output.fasta}
         gzip -dc {output.gtf}.gz > {output.gtf}.tmp && mv {output.gtf}.tmp {output.gtf}
         rm -f {output.fasta}.gz {output.gtf}.gz) > {log} 2>&1
        """


rule index_fasta:
    input:
        GENOME,
    output:
        f"{GENOME}.fai",
    shell:
        "samtools faidx {input}"


rule bwa_index:
    input:
        GENOME,
    output:
        multiext(GENOME, ".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac"),
    log:
        f"{REF_DIR}/logs/bwa_index.log",
    resources:
        mem_mb=40000,
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


rule snpeff_build:
    input:
        fasta=GENOME,
        gtf=GTF,
        codon_tables="workflow/resources/snpeff_codon_tables.config",
    output:
        db=SNPEFF_DB,
        config=SNPEFF_CONFIG,
    params:
        data_dir=lambda w: os.path.abspath(SNPEFF_DIR),
        db_dir=f"{SNPEFF_DIR}/{REF_NAME}",
        name=REF_NAME,
    log:
        f"{REF_DIR}/logs/snpeff_build.log",
    shell:
        """
        mkdir -p {params.db_dir}
        ln -sf $(readlink -f {input.fasta}) {params.db_dir}/sequences.fa
        ln -sf $(readlink -f {input.gtf}) {params.db_dir}/genes.gtf
        printf '%s\n' \
            'data.dir = {params.data_dir}' \
            '{params.name}.genome : Danio rerio' \
            '{params.name}.MT.codonTable : Vertebrate_Mitochondrial' \
            > {output.config}
        cat {input.codon_tables} >> {output.config}
        snpEff -Xmx16g build -gtf22 -noCheckCds -noCheckProtein \
            -c {output.config} -v {params.name} > {log} 2>&1
        """


rule genome_chunks:
    """Split the genome into similarly sized groups of chromosomes/scaffolds."""
    input:
        f"{GENOME}.fai",
    output:
        expand(f"{REF_DIR}/chunks/chunk_{{chunk}}.bed", chunk=CHUNKS),
    params:
        n=N_CHUNKS,
        outdir=f"{REF_DIR}/chunks",
    script:
        "../scripts/genome_chunks.py"


rule known_variants:
    """Ensembl/dbSNP zebrafish variants, split and left-aligned like our own calls.

    The Ensembl file has no ##contig lines, so they are added from the genome index.
    """
    input:
        fasta=GENOME,
        fai=f"{GENOME}.fai",
    output:
        vcf=KNOWN_VCF,
        tbi=f"{KNOWN_VCF}.tbi",
    params:
        url=config.get("known_variants", {}).get("url", ""),
    threads: 4
    log:
        f"{REF_DIR}/logs/known_variants.log",
    shell:
        """
        (wget -q -c -O {output.vcf}.download.gz {params.url}
         bcftools reheader --fai {input.fai} {output.vcf}.download.gz \
           | bcftools annotate -x INFO -Ou \
           | bcftools norm --threads {threads} -f {input.fasta} -m -any --check-ref w -Ou \
           | bcftools norm --threads {threads} -d exact -Oz -o {output.vcf}
         tabix -f -p vcf {output.vcf}
         rm -f {output.vcf}.download.gz) > {log} 2>&1
        """


rule protein_references:
    """Proteins used to check missense alleles: zebrafish (Ensembl + RefSeq) and orthologs."""
    output:
        **PROTEIN_REFS,
    params:
        ensembl_ftp=config.get("protein_check", {}).get("ensembl_ftp", ""),
        release=lambda w: config["protein_check"]["ensembl_ftp"].rstrip("/").split("release-")[-1],
        species=config.get("protein_check", {}).get("species", []),
        refseq_xref_url=config.get("protein_check", {}).get("refseq_xref_url", ""),
        refseq_pep_url=config.get("protein_check", {}).get("refseq_protein_url", ""),
    log:
        f"{REF_DIR}/logs/protein_references.log",
    script:
        "../scripts/prepare_protein_refs.py"
