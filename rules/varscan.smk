rule samtools_mpileup:
    input:
        # single or list of bam files
        bam="bam/{sample}.bam",
        reference_genome="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"
    output:
        temp("mpileup/{sample}.mpileup.gz")
    log:
        "logs/varscan_calling/mpileup/{sample}.log"
    threads:
        2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    params:
        extra="-d 10000",  # optional
    wrapper:
        "0.67.0/bio/samtools/mpileup"

rule unzip_mpileup:
    input:
        pileup_gz="mpileup/{sample}.mpileup.gz"
    output:
        "mpileup/{sample}.mpileup"
    log:
        "logs/varscan_calling/mpileup/gunzip/{sample}.log"
    threads:
        2
    shell:
        "gunzip -c -d {input} > {output} 2> {log}"



rule varscan_snp:
    input:
        "mpileup/{sample}.mpileup"
    output:
        "vc/varscan_snp/{sample}_tmp.vcf"
    message:
        "Calling SNP with Varscan2"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    threads:  # Varscan does not take any threading information
        1     # However, mpileup might have to be unzipped.
              # Keep threading value to one for unzipped mpileup input
              # Set it to two for zipped mipileup files
    log:
        "logs/varscan_calling/snp/calling/{sample}.log"
    shell:
       "java -jar ~/anaconda3/envs/variant_calling/share/varscan-2.4.4-0/VarScan.jar mpileup2snp {input} --output-vcf 1 > {output} 2> {log} "


"""
Création d'un fichier texte d'une ligne avec le nom que l'on souhait donner au "sample". Utile à la fin quand on merge les fichiers
des différents Variant Callers utilisés, afin de répere quels variants sont trouvés par quels outils
"""
rule add_varscansnp:
    output:
        temp("vc/varscan_snp/{sample}_varscan_snp.txt")  # either .vcf or .bcf
    log:
        "logs/varscan_calling/snp/rename_sample/{sample}.log"
    threads: 1
    resources:
        time_min = (
            lambda wildcars, attempt: min(10 * attempt, 20)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(128 * attempt, 512)
        )
    shell:
        "echo {wildcards.sample}_VSNP > {output}"

"""
This rule calls variants change sample name (utile pour le merge)
eg : DOOOQCL --> DOOQCL_VSNP
"""

rule change_sample_VSNP:
    input:
        vcf="vc/varscan_snp/{sample}_tmp.vcf",
        sample="vc/varscan_snp/{sample}_varscan_snp.txt",
    output:
        temp("vc/varscan_snp/{sample}.vcf")  # either .vcf or .bcf
    log:
        "logs/varscan_calling/snp/reheader/{sample}.log"
    threads: 2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    shell:
        "bcftools reheader -s {input.sample} {input.vcf} > {output} 2> {log}"

"""
Thir rule compress vcf : varscan snp
"""
rule compress_vcf_varscan_snp:
    input:
        "vc/varscan_snp/{sample}.vcf"
    output:
        "vc/varscan_snp/{sample}.vcf.gz"
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
        "logs/varscan_calling/snp/zip/{sample}.log"
    shell:
        "bgzip {input} -c --threads {threads} > {output} 2> {log}"


rule tabix_var_snp:
    input:
        "vc/varscan_snp/{sample}.vcf.gz"
    output:
        "vc/varscan_snp/{sample}.vcf.gz.tbi"
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
        "logs/varscan_calling/snp/tabix/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"


rule varscan_indel:
    input:
        "mpileup/{sample}.mpileup"
    output:
        temp("vc/varscan_indel/{sample}_tmp.vcf")
    message:
        "Calling Indel with Varscan2"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/varscan_calling/indel/calling/{sample}.log"
    shell:
        "java -jar ~/anaconda3/envs/variant_calling/share/varscan-2.4.4-0/VarScan.jar mpileup2indel {input} --output-vcf 1 > {output} 2> {log}"

"""
This rule create the file to rename the sample in the vc file for the caller freebayes
"""
rule add_varscan_indel:
    output:
        temp("vc/varscan_indel/{sample}_varscan_indel.txt")  # either .vcf or .bcf
    log:
        "logs/varscan_calling/indel/rename_sample/{sample}.log"
    threads: 1
    resources:
        time_min = (
            lambda wildcars, attempt: min(10 * attempt, 20)
        ),
        mem_mb = (
            lambda wildcars, attempt: min(128 * attempt, 512)
        )
    shell:
        "echo {wildcards.sample}_VIND > {output} 2> {log}"

"""
This rule calls variants change sample name (utile pour le merge)
eg : DOOOQCL --> DOOQCL_VSNP
"""

rule change_sample_VIND:
    input:
        vcf="vc/varscan_indel/{sample}_tmp.vcf",
        sample="vc/varscan_indel/{sample}_varscan_indel.txt",
    output:
        temp("vc/varscan_indel/{sample}.vcf")  # either .vcf or .bcf
    log:
        "logs/varscan_calling/indel/reheader/{sample}.log"
    threads: 2
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    shell:
        "bcftools reheader -s {input.sample} {input.vcf} > {output} 2> {log}"


rule compress_vcf_varscan_indel:
    input:
        "vc/varscan_indel/{sample}.vcf"
    output:
        "vc/varscan_indel/{sample}.vcf.gz"
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
        "logs/varscan_calling/indel/zip/{sample}.log"
    shell:
        "bgzip {input} -c --threads {threads} > {output} 2> {log}"


rule tabix_varscan_indel:
    input:
        "vc/varscan_indel/{sample}.vcf.gz"
    output:
        "vc/varscan_indel/{sample}.vcf.gz.tbi"
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
        "logs/varscan_calling/indel/tabix/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"
