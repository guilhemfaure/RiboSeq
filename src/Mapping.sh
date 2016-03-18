rRNA=../../../assembly/Bsubtilis/NC_00964/tRNA_rRNA # index of tRNA_rRNA sequence
genome=../../../assembly/Bsubtilis/NC_00964/genome # index of the genome sequence
read=SRR407278_clip.fastq.gz # read free from adapter

# Remove tRNA and rRNA from reads
hisat2 \
-x $rRNA \
-q ${read%.*} \
--phred33 \
--un noRNA.fastq \
-p 10 \
-t \
 > /dev/null


read2=../noRNA.fastq
sam=NC_00964.sam
# Map reads on genome
hisat2 \
-x ../$genome \
-q $read2 \
--phred33 \
-p 10 \
-t \
-S $sam


samtools view -bS NC_000913.2.sam > NC_000913.2.bam
samtools sort NC_000913.2.bam NC_000913.2.sort
samtools index NC_000913.2.sort.bam


# Then see RibosomeOccupancy.py

