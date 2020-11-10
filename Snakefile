import glob

SAMPLES=[]

sample_id_list = glob.glob('/home/claire/Documents/bam/*.bam')
for name in sample_id_list:
    SAMPLES.append(name.split('/')[5].split('.')[0])
print(SAMPLES)

genome="/home/claire/Documents/genome/Homo_sapiens_GRCh38.dna.primary_assembly.fa"

include: "rules/index.smk"
include: "rules/freebayes.smk"
include: "rules/varscan.smk"
include: "rules/haplotype_caller.smk"
include: "rules/merge.smk"
include: "rules/deepvariant.smk"
include: "rules/copy.smk"


rule all:
    input:
        test=expand(
            "vc/hc/{sample}_tmp.vcf",
            sample=SAMPLES
        )
    message:
        "Finish"
