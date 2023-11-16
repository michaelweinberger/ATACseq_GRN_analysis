#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=1-00:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=12

module load bowtie/1.2.3
module load samtools/1.10
module load bedtools/2.29.2
module load ucsctools/385
module load python-cbrg


#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################



############################################################



### concatenate fastq files across different sequencing lanes
cat_fastq_dir="${out_dir}/cat_fastq_ATAC"
[ ! -d "$cat_fastq_dir" ] && mkdir -p "$cat_fastq_dir"

# loop through list of fastq directories
len=${#ATAC_fastq_dirs[@]}

for (( k=0; k<$len; k++ ))
do
	cd ${ATAC_fastq_dirs[k]}

	# loop over fastq files in directory
	for i in $(ls *.fastq.gz | rev | cut -c 22- | rev | uniq)
	do
		echo Concatenating ${conditions[k]} ${i} fastq files

		cat ${i}*R1* > ${cat_fastq_dir}/${conditions[k]}_${i}_all_lanes_R1.fastq.gz
		cat ${i}*R2* > ${cat_fastq_dir}/${conditions[k]}_${i}_all_lanes_R2.fastq.gz
		wait
	done
done

wait




### run fastqc quality check
#fastqc_dir="${out_dir}/cat_fastq_fastqc"
#[ ! -d "$fastqc_dir" ] && mkdir -p "$fastqc_dir"
#cd $cat_fastq_dir
#fastqc *
#wait
#mv *fastqc* $fastqc_dir




### map reads
sam_dir="${out_dir}/bowtie_${genome}_mapped"
[ ! -d "$sam_dir" ] && mkdir -p "$sam_dir"

cd $cat_fastq_dir


# replace dashes and dots in concatenated fastq file names with underscores
for i in $(ls *.fastq.gz)
do
	new_name=$(echo ${i} | rev | cut -c 10- | rev | uniq | tr '-' '_' | tr '.' '_')
	mv ${i} ${new_name}.fastq.gz

done


# loop over concatenated fastq files
for i in $(ls *.fastq.gz | rev | cut -c 23- | rev | uniq)
do
	if [ ! -f "${sam_dir}/${i}_bowtie.log" ]; then
		echo Mapping ${i}
		gunzip ${i}_all_lanes_R1.fastq.gz
		gunzip ${i}_all_lanes_R2.fastq.gz
		wait
		bowtie --sam --maxins 1250 -p 12  ${genome_dir}/${genome} -1 ${i}_all_lanes_R1.fastq -2 ${i}_all_lanes_R2.fastq 1> ${sam_dir}/${i}.sam 2> ${sam_dir}/${i}_bowtie.log
		wait
		gzip ${i}_all_lanes_R1.fastq
		gzip ${i}_all_lanes_R2.fastq
		wait
	fi
done

wait




### process sam files
bam_dir="${sam_dir}/BAM_files"
[ ! -d "$bam_dir" ] && mkdir -p "$bam_dir"

cd $sam_dir

for j in $(ls *.sam | rev | cut -c 5- | rev | uniq)
do
	# only process SAM file if final BAM file is not already present
	if [ ! -f "${bam_dir}/${j}.rmdup.bam" ]; then

		# compress sam files into bam format
		echo Converting ${j} SAM file to BAM file
		samtools view -bS -o ${bam_dir}/${j}.bam ${j}.sam 

		# sort bam files
		echo Sorting  ${j} BAM file
		samtools sort -o ${bam_dir}/${j}.sorted.bam ${bam_dir}/${j}.bam

		# remove duplicate reads from bam files
		echo Removing duplicates from ${j} sorted BAM file
		samtools rmdup ${bam_dir}/${j}.sorted.bam ${bam_dir}/${j}.rmdup.bam

		# index bam files
		echo Indexing ${j} BAM file
		samtools index ${bam_dir}/${j}.rmdup.bam

		# remove unnecessary files
		rm ${j}.sam
		rm ${bam_dir}/${j}.bam
		rm ${bam_dir}/${j}.sorted.bam
	fi
done

wait




### count entries in final bam file
cd $bam_dir

for j in $(ls *.rmdup.bam | rev | cut -c 11- | rev | uniq)
do
	echo Number of reads in $j sorted + no duplicates BAM file:
	samtools view -c ${j}.rmdup.bam
done

wait




### generate bigWig files for UCSC genome browser session
for j in $(ls *.rmdup.bam | rev | cut -c 11- | rev | uniq)
do
	# convert bam file into bedgraph file
	if [ ! -f "${bam_dir}/${j}.bg" ]; then
		echo Converting ${j} BAM file to BEDGRAPH file
		genomeCoverageBed -bg -split -ibam ${bam_dir}/${j}.rmdup.bam > ${bam_dir}/${j}.bg 
		bedSort ${bam_dir}/${j}.bg ${bam_dir}/${j}.bg
	fi
	
	# convert bedgraph file into bigwig file
	if [ ! -f "${bam_dir}/${j}.bw" ]; then
		echo Converting ${j} BEDGRAPH file to BIGWIG file
		bedGraphToBigWig ${bam_dir}/${j}.bg ${genome_dir}/${genome}.chrom.sizes ${bam_dir}/${j}.bw
	fi
done

wait



echo All done! 

