"""
This rule creates a fasta sequence index.
More information: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/samtools/faidx.html
"""
rule samtools_faidx:
    input:
        "genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"
    output:
        "genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa.fai"
    message:
        "Building fasta-sequence index for reference genome sequence"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 20, 200)
        )
    log:
        "logs/index/faidx.log"
    wrapper:
        "0.67.0/bio/samtools/faidx"



"""
This rule creates a fasta sequence dictionnary.
More information: https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/picard/createsequencedictionary.html
"""
rule create_dict:
    input:
        "genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"
    output:
        "genome/Homo_sapiens_GRCh38.dna.primary_assembly.dict"
    message:
        "Building fasta-sequence dictionnary for reference genome sequence"
    threads:
        1
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 1024, 10240)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 10, 200)
        )
    log:
        "logs/index/createseqdico.log"
    wrapper:
        "0.67.0/bio/picard/createsequencedictionary"
