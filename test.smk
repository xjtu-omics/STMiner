STAR = "/data/home/pssun/miniconda2/envs/ngs/bin/STAR"
dir_work = "/data/home/pssun/data/Zhang/"
INDEX = ["D3-1", "D3-2", "D3-3"]
rule all:
    input:  expand(dir_work + "bulk_bam_data/D3/{index}_Aligned.sortedByCoord.out.bam", index=INDEX)

rule map_to_mum:
    input:
        r1 = dir_work + "X101SC23011015-Z02-J001/00.CleanData/{index}/{index}_1.clean.fq.gz", 
        r2 = dir_work + "X101SC23011015-Z02-J001/00.CleanData/{index}/{index}_2.clean.fq.gz"
    output: 
        dir_work + "bulk_bam_data/D3/{index}_Aligned.sortedByCoord.out.bam"
    shell:
        "{STAR} --runThreadN 24 \
        --runMode alignReads \
        --readFilesCommand zcat \
        --quantMode TranscriptomeSAM GeneCounts \
        --twopassMode Basic \
        --outSAMtype BAM SortedByCoordinate \
        --genomeDir ~/data/Zhang/genome/mum \
        --readFilesIn  {input} \
        --outFileNamePrefix {dir_work}bulk_bam_data/D3/{index}/{wildcards.index}_"
        