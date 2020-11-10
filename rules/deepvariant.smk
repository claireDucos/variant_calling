"""
Variant Calling using DeepVariant
https://github.com/google/deepvariant
"""
rule deepvariant_calling:
    input:
        bam="bam/{sample}.bam",
        ref="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa",
        index="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa.fai"
    output:
        vcf=temp("vc/deep/{sample}_tmp.vcf")
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    params:
        model="wes",   # {wgs, wes}
        extra=""
    threads: 2
    log:
        "logs/deepvariant/{sample}.log"
    wrapper:
        "0.67.0/bio/deepvariant"

"""
Création d'un fichier texte d'une ligne avec le nom que l'on souhait donner au "sample". Utile à la fin quand on merge les fichiers
des différents Variant Callers utilisés, afin de répere quels variants sont trouvés par quels outils
"""
rule create_sample_change_file_deep:
    output:
        temp("vc/deep/{sample}_DEEP.txt")  # either .vcf or .bcf
    log:
        "logs/deepvariant_calling/rename_sample/{sample}_txt.log"
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512 + 512, 2048)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 120)
        )
    shell:
        "echo {wildcards.sample}_DEEP > {output}"

"""
This rule calls variants change sample name (utile pour le merge)
eg : DOOOQCL --> DOOQCL_DEEP
"""

rule rename_sample_deep:
    input:
        vcf="vc/deep/{sample}_tmp.vcf",
        sample="vc/deep/{sample}_DEEP.txt",
    output:
        temp("vc/deep/{sample}.vcf") # either .vcf or .bcf
    log:
        "logs/deepvariant_calling/reheader/{sample}.log"
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
This rule compress the vcf file : Deepvariant
"""
rule compress_vcf_deep:
    input:
        "vc/deep/{sample}.vcf"
    output:
        "vc/deep/{sample}.vcf.gz"
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
        "logs/deepvariant_calling/zip/{sample}.log"
    shell:
        "bgzip {input} -c --threads {threads} > {output}"


"""
This rule index the vcf file : Deepvariant
"""
rule tabix_deep:
    input:
        "vc/deep/{sample}.vcf.gz"
    output:
        "vc/deep/{sample}.vcf.gz.tbi"
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
        "logs/deepvariant_calling/tabix/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"
