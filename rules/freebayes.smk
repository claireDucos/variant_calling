"""
This rule calls variants with FreeBayes
More information at:
https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/bio/freebayes.html
"""

rule freebayes_variant_calling:
    input:
        samples="bam/{sample}.bam",
        indexes="bam/{sample}.bam.bai",
        ref="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa",
        index="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa.fai"

    output:
        temp("vc/freebayes/{sample}_tmp.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes_calling/calling/{sample}.log"
    params:
        extra="",         # optional parameters
        chunksize=100000  # reference genome chunk size for parallelization (default: 100000)
    threads: 2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 128, 512)
        )
    wrapper:
        "0.67.0/bio/freebayes"

"""
Création d'un fichier texte d'une ligne avec le nom que l'on souhait donner au "sample". Utile à la fin quand on merge les fichiers
des différents Variant Callers utilisés, afin de répere quels variants sont trouvés par quels outils
"""
rule create_sample_change_file_FB:
    output:
        temp("vc/freebayes/{sample}_freebayes.txt")  # either .vcf or .bcf
    log:
        "logs/freebayes_calling/rename_sample/{sample}_txt.log"
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512 + 512, 2048)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 120)
        )
    shell:
        "echo {wildcards.sample}_FB > {output}"

"""
This rule calls variants change sample name (utile pour le merge)
eg : DOOOQCL --> DOOQCL_FB
"""

rule rename_sample_FB:
    input:
        vcf="vc/freebayes/{sample}_tmp.vcf",
        sample="vc/freebayes/{sample}_freebayes.txt",
    output:
        temp("vc/freebayes/{sample}.vcf")  # either .vcf or .bcf
    log:
        "logs/freebayes_calling/reheader/{sample}.log"
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512 + 512, 2048)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 120)
        )
    shell:
        "bcftools reheader -s {input.sample} {input.vcf} > {output} 2> {log}"

"""
This rule compress vcf file : Freebayes
"""
rule compress_vcf_FB:
    input:
        "vc/freebayes/{sample}.vcf"
    output:
        "vc/freebayes/{sample}.vcf.gz"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    threads:
        2
    log:
        "logs/freebayes_calling/zip/{sample}.log"
    shell:
        "bgzip {input} -c --threads {threads} > {output}"

"""
This rule index the vcf file : FreeBayes
"""
rule tabix_fb:
    input:
        "vc/freebayes/{sample}.vcf.gz"
    output:
        "vc/freebayes/{sample}.vcf.gz.tbi"
    params:
        # pass arguments to tabix (e.g. index a vcf)
        "-p vcf"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/freebayes_calling/tabix/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"
