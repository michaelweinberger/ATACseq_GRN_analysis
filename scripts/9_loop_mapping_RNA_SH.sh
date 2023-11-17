#!/bin/sh

# Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-12:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12



#########################################################
### This script:
# creates STAR genome index
# performs STAR mapping of bulk RNA-seq data
# creates BAM files, sorts and removes duplicate reads
# counts number of entries in BAM files
# performs featureCounts
#########################################################



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

#################################################################################################################################################


#################################################################################################################################################



module load STAR/2.7.10b
module load samtools/1.10
module load bedtools/2.29.2
module load ucsctools/385
module load subread/2.0.6



### generate genome files
# check if directory for genome files exists, if not create it using mkdir
[ ! -d "$genome_dir" ] && mkdir -p "$genome_dir"


# download fasta file if not present
if [ ! -f "${genome_dir}/${genome}.fa" ]; then
	rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/${species_latin}/dna/${species_latin_2}.${genome}.dna.toplevel.fa.gz $genome_dir
	rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/release-110/fasta/${species_latin}/dna_index/${species_latin_2}.${genome}.dna.toplevel.fa.gz.fai $genome_dir

	#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/current_fasta/${species_latin}/dna/${species_latin_2}.${genome}.dna.toplevel.fa.gz $genome_dir
	#rsync -avzP rsync://ftp.ensembl.org/ensembl/pub/current_fasta/${species_latin}/dna_index/${species_latin_2}.${genome}.dna.toplevel.fa.gz.fai $genome_dir

	mv ${genome_dir}/${species_latin_2}.${genome}.dna.toplevel.fa.gz ${genome_dir}/${genome}.fa.gz
	mv ${genome_dir}/${species_latin_2}.${genome}.dna.toplevel.fa.gz.fai ${genome_dir}/${genome}.fa.fai
	gunzip ${genome_dir}/${genome}.fa.gz
fi


# generate chromosome sizes file if not present
if [ ! -f "${genome_dir}/${genome}.chrom.sizes" ]; then
	awk -F '\t' '{print $1, $2}' ${genome_dir}/${genome}.fa.fai > ${genome_dir}/${genome}.chrom.sizes
fi


# download gtf file if not present
if [ ! -f "${genome_dir}/${genome}.gtf" ]; then
	rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-110/gtf/${species_latin}/${species_latin_2}.${genome}.110.gtf.gz $genome_dir
	mv ${genome_dir}/${species_latin_2}.${genome}.110.gtf.gz ${genome_dir}/${genome}.gtf.gz
	gunzip ${genome_dir}/${genome}.gtf.gz
fi



# subset genome files to main chromosomes, as there might be many additional contigs present, meaning large memory usage in downstream analysis 
# subset chromosome sizes file to main chromosomes
if [ ! -f "${genome_dir}/${genome}_subset.chromosomes.txt" ]; then
	head -n $chrom_number ${genome_dir}/${genome}.chrom.sizes > ${genome_dir}/${genome}_subset.chrom.sizes
	awk 'BEGIN { OFS = "\n" }{print $1}' ${genome_dir}/${genome}_subset.chrom.sizes > ${genome_dir}/${genome}_subset.chromosomes.txt
fi

# subset genome fasta file to main chromosomes
if [ ! -f "${genome_dir}/${genome}_subset.fa" ]; then
	samtools faidx ${genome_dir}/${genome}.fa -r ${genome_dir}/${genome}_subset.chromosomes.txt > ${genome_dir}/${genome}_subset.fa
	samtools faidx ${genome_dir}/${genome}_subset.fa > ${genome_dir}/${genome}_subset.fa.fai
fi

# subset gtf file to main chromosomes
if [ ! -f "${genome_dir}/${genome}_subset.gtf" ]; then

	head -n 5 ${genome_dir}/${genome}.gtf > ${genome_dir}/${genome}_subset.gtf

	for k in $(awk 'BEGIN { OFS = "\n" }{print $1}' ${genome_dir}/${genome}_subset.chromosomes.txt)
	do
		awk -F $'\t' -v k="${k}" '$1==k' ${genome_dir}/${genome}.gtf > ${genome_dir}/${genome}_subset_${k}.gtf
		#wc -l ${genome_dir}/${genome}_subset_${k}.gtf
	done
	cat ${genome_dir}/${genome}_subset_*.gtf >> ${genome_dir}/${genome}_subset.gtf
	#wc -l ${genome_dir}/${genome}_subset.gtf

	rm ${genome_dir}/${genome}_subset_*.gtf
fi
	
genome_size=$(awk '{SUM+=$2}; END {print SUM}' ${genome_dir}/${genome}_subset.chrom.sizes)



wait




### create genome index for STAR if not present
# check if directory for genome index exists, if not create it using mkdir
star_dir="${genome_dir}/STAR_${genome}_subset"
[ ! -d "$star_dir" ] && mkdir -p "$star_dir"

if [ ! -f "${star_dir}/SA" ]; then
	STAR --runThreadN 12 \
    	--runMode genomeGenerate \
   	--genomeDir ${star_dir} \
    	--genomeFastaFiles ${genome_dir}/${genome}_subset.fa \
   	--sjdbGTFfile ${genome_dir}/${genome}_subset.gtf
fi

wait



### concatenate fastq files across different sequencing lanes
cat_fastq_dir_2="${out_dir}/cat_fastq_RNA_ananse"
[ ! -d "$cat_fastq_dir_2" ] && mkdir -p "$cat_fastq_dir_2"

# loop through list of fastq directories
len=${#RNA_fastq_dirs_ananse[@]}

for (( k=0; k<${len}; k++ ))
do
	cd ${RNA_fastq_dirs_ananse[k]}

	# loop over fastq files in directory
	for i in $(ls *.fastq.gz | rev | cut -c 22- | rev | uniq)
	do
		echo Concatenating ${conditions_ananse_RNA[k]} ${i} fastq files

		cat ${i}*R1* > ${cat_fastq_dir_2}/${conditions_ananse_RNA[k]}_${i}_all_lanes_R1.fastq.gz
		cat ${i}*R2* > ${cat_fastq_dir_2}/${conditions_ananse_RNA[k]}_${i}_all_lanes_R2.fastq.gz
		wait
	done
done

wait




### map reads and process mapped files

sam_dir_2="${out_dir}/STAR_${genome}_mapped_ananse"
[ ! -d "$sam_dir_2" ] && mkdir -p "$sam_dir_2"

cd $cat_fastq_dir_2


# replace dashes and dots in concatenated fastq file names with underscores
for i in $(ls *.fastq.gz)
do
	new_name=$(echo ${i} | rev | cut -c 10- | rev | uniq | tr '-' '_' | tr '.' '_')
	mv ${i} ${new_name}.fastq.gz

done


# loop over concatenated fastq files
for i in $(ls *.fastq.gz | rev | cut -c 23- | rev | uniq)
do

	# only do mapping if log file does not already exist
	if [ ! -f "${sam_dir_2}/${i}_mapped/${i}_Log.final.out" ]; then

		mkdir ${sam_dir_2}/${i}_mapped/
		echo Mapping $i
		STAR --runThreadN 12 \
  		--genomeDir ${star_dir} \
  		--readFilesIn ${i}_all_lanes_R1.fastq.gz ${i}_all_lanes_R2.fastq.gz \
  		--readFilesCommand zcat \
  		--outFileNamePrefix ${sam_dir_2}/${i}_mapped/${i}_ \
  		--outSJfilterCountUniqueMin -1 -1 -1 -1 \
  		--outSJfilterCountTotalMin 30 10 10 10 \
  		--outSAMstrandField intronMotif \
  		--outFilterIntronMotifs RemoveNoncanonical
	fi
done

wait




### process sam files
bam_dir_2="${sam_dir_2}/BAM_files"
[ ! -d "$bam_dir_2" ] && mkdir -p "$bam_dir_2"

for j in $(find $sam_dir_2 -name '*_Aligned.out.sam' | awk -F/ '{print $(NF-0)}' | rev | cut -c 17- | rev | uniq)
do

	# only process SAM file if final BAM file is not already present
	if [ ! -f "${bam_dir_2}/${j}.rmdup.bam" ]; then

		# compress sam files into bam format
		echo Converting ${j} SAM file to BAM file
		samtools view -bS -o ${bam_dir_2}/${j}.bam ${sam_dir_2}/${j}_mapped/${j}_Aligned.out.sam 

		# sort bam files
		echo Sorting  ${j} BAM file
		samtools sort -o ${bam_dir_2}/${j}.sorted.bam ${bam_dir_2}/${j}.bam

		# remove duplicate reads from bam files
		echo Removing duplicates from ${j} sorted BAM file
		samtools rmdup ${bam_dir_2}/${j}.sorted.bam ${bam_dir_2}/${j}.rmdup.bam

		# remove unnecessary files
		rm ${sam_dir_2}/${j}_mapped/${j}_Aligned.out.sam
		rm ${bam_dir_2}/${j}.bam
		rm ${bam_dir_2}/${j}.sorted.bam
	fi
done




### count entries in final bam file
for j in $(ls ${bam_dir_2}/*.rmdup.bam | uniq)
do
	echo Number of reads in ${j}:
	samtools view -c ${j}
done

wait




### generate bigWig files for UCSC genome browser session
for j in $(ls ${bam_dir_2}/*.rmdup.bam | rev | cut -c 11- | rev | uniq)
do
	# convert bam file into bedgraph file
	if [ ! -f "${j}.bg" ]; then
		echo Converting ${j} BAM file to BEDGRAPH file
		genomeCoverageBed -bg -split -ibam ${j}.rmdup.bam > ${j}.bg 
		bedSort ${j}.bg ${j}.bg
	fi
	
	# convert bedgraph file into bigwig file
	if [ ! -f "${j}.bw" ]; then
		echo Converting ${j} BEDGRAPH file to BIGWIG file
		bedGraphToBigWig ${j}.bg ${genome_dir}/${genome}_subset.chrom.sizes ${j}.bw 
	fi
done

wait




### generate gene expression count matrix
feature_input=$(ls ${bam_dir_2}/*.rmdup.bam)

if [ ! -f "${sam_dir_2}/featureCounts_raw.txt" ]; then
	featureCounts -p -B -a ${genome_dir}/${genome}_subset.gtf \
	-t exon \
	-g gene_id \
	-o ${sam_dir_2}/featureCounts_raw.txt \
	$feature_input
fi



echo All done!
