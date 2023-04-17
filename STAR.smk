STAR = "/data/home/pssun/miniconda2/envs/ngs/bin/STAR"
dir_work = "/data/home/pssun/data/Zhang/"
INDEX = ["ctrl-1","ctrl-2","ctrl-3","D3-1","D3-2","D3-3","D1-4","D1-2","D1-3","D5-1","D5-2","D5-7","D7-1","D7-2","D7-3"]
ITEM = [ dir_work + f"new_bulk_bam_data/{i.split('-')[0]}/{i}_Aligned.sortedByCoord.out.bam" for i in INDEX]
rule all:
    input: ITEM
rule map_to_mum:
    input:
        r1 = dir_work + "new_data/00.CleanData/{day}-{sample}/{day}-{sample}_1.clean.fq.gz",
        r2 = dir_work + "new_data/00.CleanData/{day}-{sample}/{day}-{sample}_2.clean.fq.gz"
    output:
        dir_work + "new_bulk_bam_data/{day}/{day}-{sample}_Aligned.sortedByCoord.out.bam"
    shell:
        "{STAR} --runThreadN 24 \
        --runMode alignReads \
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir ~/data/Zhang/genome/mum \
        --readFilesIn {input} \
        --outFileNamePrefix {dir_work}new_bulk_bam_data/{wildcards.day}/{wildcards.day}-{wildcards.sample}/mum_{wildcards.day}-{wildcards.sample}_"