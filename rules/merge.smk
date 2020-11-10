"""
This rule merge the vcf frome the different variant callers
"""
rule bcftools_merge:
    input:
        freebayes="vc/freebayes/{sample}.vcf.gz",
        hc="vc/hc/{sample}.vcf.gz",
        vsnp="vc/varscan_snp/{sample}.vcf.gz",
        vindel="vc/varscan_indel/{sample}.vcf.gz",
        deep="vc/deep/{sample}.vcf.gz",
        index_hc="vc/hc/{sample}.vcf.gz.tbi",
        index_freebayes="vc/freebayes/{sample}.vcf.gz.tbi",
        index_vsnp="vc/varscan_snp/{sample}.vcf.gz.tbi",
        index_vindel="vc/varscan_indel/{sample}.vcf.gz.tbi",
        index_deep="vc/deep/{sample}.vcf.gz.tbi"
    output:
        "vcf_merge/{sample}.vcf"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 200)
        )
    threads:
        2
    log:
        "logs/vcf_merge/{sample}.vcf"
    shell:
        "bcftools merge  {input.freebayes} {input.hc} {input.vsnp} {input.vindel} {input.deep} -o {output} --force-sample 2> {log}"

"""
This rule split multiallelic site in the merged vcf
"""
rule bcftools_split_multiallelic_site:
    input:
        vcf="vcf_merge/{sample}.vcf"
    output:
        multi=temp("vcf_split_multi/{sample}.vcf")
    message:
        "Split multi allelic site ie A,AT,TAT"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/bcftools/split/{sample}.log"
    threads:
        2
    shell:
        "bcftools norm -m-any {input} > {output} 2> {log}"


"""
This rule normalise variants in the vcf file
"""

rule bcftools_normalize_variant:
    input:
        vcf="vcf_split_multi/{sample}.vcf",
        ref="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"
    output:
        norm="final_vcf/{sample}.vcf.gz"
    message:
        "normalize variant"
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/bcftools/normalize/{sample}.log"
    threads:
        2
    shell:
        "bcftools norm -f {input.ref} {input.vcf} -O z > {output} 2> {log}"

"""
This rule index the vcf file : merged vcf
"""
rule tabix_final:
    input:
        "final_vcf/{sample}.vcf.gz"
    output:
        "final_vcf/{sample}.vcf.gz.tbi"
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
        "logs/tabix_merge/{sample}.log"
    wrapper:
        "0.67.0/bio/tabix"
