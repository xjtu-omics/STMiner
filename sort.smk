SAMTOOLS = "/data/home/pssun/miniconda2/envs/ngs/bin/samtools"
dir_work = "/data/home/pssun/data/Zhang/"
INDEX = ["ctrl-1","ctrl-2","ctrl-3","D3-1","D3-2","D3-3","D1-4","D1-2","D1-3","D5-1","D5-2","D5-7","D7-1","D7-2","D7-3"]

ITEM = [dir_work + f"new_bulk_bam_data/sorted/mum_{i}_Aligned.sortedByCoord.out.sorted.bam" for i in INDEX]

rule all:
    input: ITEM
rule build_index:
    input:
        r = dir_work + "new_bulk_bam_data/{day}/{day}-{sample}/mum_{day}-{sample}_Aligned.sortedByCoord.out.bam"
    output:
        dir_work + "new_bulk_bam_data/sorted/mum_{day}-{sample}_Aligned.sortedByCoord.out.sorted.bam"
    shell:
        "{SAMTOOLS} sort \
        -n {input} \
        -o {dir_work}new_bulk_bam_data/sorted/mum_{wildcards.sample}_Aligned.sortedByCoord.out.sorted.bam"