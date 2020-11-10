"""
Variant Caller using Haplotype Caller from Gatk
https://gatk.broadinstitute.org/hc/en-us/articles/360035531412-HaplotypeCaller-in-a-nutshell
"""
rule haplotypeCaller_calling:
    input:
        ref = "genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa",
        bam = "bam/{sample}.bam",
        dic="genome/Homo_sapiens_GRCh38.dna.primary_assembly.dict",
        fai="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa.fai"
    output:
        vcf = temp("vc/hc/{sample}_tmp.vcf")
    message:
        "Calling variants with Haplotype Caller on {wildcards.sample}"
    threads:
        4
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 128, 512)
        )
    log:
        "logs/haplotypeCaller_calling/calling/{sample}.log"
    shell:
        """gatk --java-options "-Xmx8G" HaplotypeCaller -R {input.ref} -I {input.bam} -O {output.vcf} --sequence-dictionary {input.dic} --create-output-variant-index FALSE 2> {log} """


"""
Création d'un fichier texte d'une ligne avec le nom que l'on souhait donner au "sample". Utile à la fin quand on merge les fichiers
des différents Variant Callers utilisés, afin de répere quels variants sont trouvés par quels outils
"""
rule create_sample_change_file_HC:
    output:
        temp("vc/hc/{sample}_hc.txt")  # either .vcf or .bcf
    log:
        "logs/haplotypeCaller_calling/rename_sample/{sample}_txt.log"
    threads: 1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 512 + 512, 2048)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 120)
        )
    shell:
        "echo {wildcards.sample}_HC > {output}"

"""
This rule calls variants change sample name (utile pour le merge)
eg : DOOOQCL --> DOOQCL_HC
"""

rule rename_sample_HC:
    input:
        vcf="vc/hc/{sample}_tmp.vcf",
        sample="vc/hc/{sample}_hc.txt",
    output:
        temp("vc/hc/{sample}.vcf") # either .vcf or .bcf
    log:
        "logs/haplotypeCaller_calling/reheader/{sample}.log"
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
This rule compress vcf file : haplotype_caller
"""
rule compress_vcf_HC:
    input:
        "vc/hc/{sample}.vcf"
    output:
        "vc/hc/{sample}.vcf.gz"
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
        "logs/haplotypeCaller_calling/zip/{sample}.log"
    shell:
        "bgzip {input} -c --threads {threads} > {output} 2> {log}"

"""
This rule index vcf file : haplotype_caller
"""
rule tabix_hc:
    input:
        "vc/hc/{sample}.vcf.gz"
    output:
        "vc/hc/{sample}.vcf.gz.tbi"
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
    threads:
        2
    log:
        "logs/haplotypeCaller_calling/tabix/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"
