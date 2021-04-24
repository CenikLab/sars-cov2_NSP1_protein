FASTQDIR=/data/projects/NSP/runs/20210401/intermediates/rnaseq/clip
BOWTIEREF=/data/projects/NSP/revision_analysis/reference/vector_reference/nsp1_nsp2_and_transcriptome
BOWTIEPARAMS="--threads 8 -L 15 "


experiments=`ls $FASTQDIR/*gz`

for e in $experiments;
do
	this_file=$e
	this_name=`basename $e | awk -F "." '{print($1)}'`
	echo $this_file
        echo $this_name
	echo "################"
        bowtie2 $BOWTIEPARAMS  -x $BOWTIEREF -q ${this_file}  2> ${this_name}.mapping.log |  samtools view -bS -q 2 -@ 8 - | samtools sort -@ 8 -o ${this_name}.bam
        samtools index -@ 8  ${this_name}.bam
        samtools idxstats -@ 8 ${this_name}.bam > ${this_name}.stats
done

