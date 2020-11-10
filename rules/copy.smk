"""
Cette règle copie de manière temporaire les fichiers bam nécessaire dans le dossier de travail
Il faut changer le chemin donc (à modifier plus tard)
"""
rule copy_bam:
    input:
        bam = "/home/claire/Documents/bam/{sample}.bam",
        bai = "/home/claire/Documents/bam/{sample}.bam.bai"
    output:
        bam = "bam/{sample}.bam",
        bai = "bam/{sample}.bam.bai"

    message:
        "Copying bam and bai files"
    threads:
        4
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/cp/{sample}.log"
    run:
        shell("cp {input.bam} {output.bam}")
        shell("cp {input.bai} {output.bai}")

"""
Cette règle copie de manière temporaire le fichier fasta du génome de ref nécessaire dans le dossier de travail
"""
rule copy_genome:
    input:
        fasta = genome
    output:
        fasta="genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"

    message:
        "Copying bam and bai files"
    threads:
        4
    resources:
        mem_mb = (
            lambda wildcards, attempt: min(attempt * 8192, 16384)
        ),
        time_min = (
            lambda wildcards, attempt: min(attempt * 45, 360)
        )
    log:
        "logs/cp/genome.log"
    shell:
        "cp {input} {output}"
