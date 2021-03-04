FASTQDIR=/data/projects/NSP/revision_analysis/external_data_download/fastqs/fastqs
BOWTIEREF=/data/projects/NSP/revision_analysis/reference/viral_reference/virus_nsp1_and_2_and_transcriptome
BOWTIEPARAMS="--threads 8 -L 15 --trim5 12 --trim3 12 "


experiments=`ls $FASTQDIR | grep "_1.fastq.gz"`

for e in $experiments;
do
	echo $e
	this_name=`echo $e | awk -F "_" '{print($1)}'`
	echo $this_name
	bowtie2 $BOWTIEPARAMS  -x $BOWTIEREF -q $FASTQDIR/$e  2> ${this_name}.mapping.log |  samtools view -bS -q 2 -@ 8 - | samtools sort -@ 8 -o ${this_name}.bam
	samtools index -@ 8  ${this_name}.bam
	samtools idxstats -@ 8 ${this_name}.bam > ${this_name}.stats
	
done
