#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-06:00:00
#SBATCH --mem=20G
#SBATCH --cpus-per-task=12



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

########################################################################################################################


########################################################################################################################



module load homer/20201202
# This command also loads requirement: R-base/4.2.0 gsl/2.6 hdf5/1.10.7 gdal/3.4.1 sqlite/3.37.2 proj/8.1.1 geos/3.9.0 jags/4.3.0 java/17.0.1 cmake/3.17.2 R-cbrg/current

module load samtools/1.10
module load ucsctools/385
module load bedtools/2.29.2




### create directory
homer_dir="${out_dir}/homer_peak_annotation"
[ ! -d "$homer_dir" ] && mkdir -p "$homer_dir"




### perform TSS enrichment analysis
# generate tag directories to use with Homer
cd $bam_dir

for j in $(ls *.rmdup.bam | rev | cut -c 11- | rev | uniq)
do
	if [ ! -d "${homer_dir}/homer_tag_dir_${j}" ]; then
		echo Generating Homer tag directory for ${j}.rmdup.bam
		makeTagDirectory ${homer_dir}/homer_tag_dir_${j} ${j}.rmdup.bam -single
	fi
done

wait


# generate histograms of tag densities across tss
tag_dirs=(${homer_dir}/homer_tag_dir_*)

echo Generating histograms of tag densities across TSS
annotatePeaks.pl tss ${genome_dir}/${genome}_subset.fa \
  -gtf ${genome_dir}/${genome}_subset.gtf \
  -d ${tag_dirs[@]} \
  -hist 10 -size 4000 > ${homer_dir}/Tag_density_tss.txt




### annotate Diffbind consensus peak set to closest genes using Homer
# annotate peaks to all genes in gtf file
echo Annotating peaks in ${diffbind_dir}/Diffbind_consensus_peak_set_for_Homer.bed to closest gene in ${genome_dir}/${genome}_subset.gtf
echo Saving as ${homer_dir}/Diffbind_consensus_peak_set_annotated_${genome}_subset.txt

annotatePeaks.pl ${diffbind_dir}/Diffbind_consensus_peak_set_for_Homer.bed \
  ${genome_dir}/${genome}_subset.fa -gtf ${genome_dir}/${genome}_subset.gtf \
  -CpG > ${homer_dir}/Diffbind_consensus_peak_set_annotated_${genome}_subset.txt


# OPTIONAL: subset gtf file to expressed genes
len=${#RNA_fastq_dirs[@]}

if [ $len > 0 ]; then
	for k in $(ls ${sam_dir_2}/Expressed_genes_*.txt | rev | cut -c 5- | rev)
	do
		# subset gtf file
		if [ ! -f "${k}.gtf" ]; then
			echo Subsetting ${genome} GTF file to genes contained in ${k}.txt
			while read line; do grep -wF "$line" ${genome_dir}/${genome}_subset.gtf; done < ${k}.txt > ${k}.gtf
		fi

		wait

		# annotate peaks to expressed genes
		gtf_name=${k##*/}
		if [ ! -f "${homer_dir}/Diffbind_consensus_peak_set_annotated_${gtf_name}.txt" ]; then
			
			echo Annotating peaks in ${diffbind_dir}/Diffbind_consensus_peak_set_for_Homer.bed to closest gene in ${k}.gtf
			echo Saving as ${homer_dir}/Diffbind_consensus_peak_set_annotated_${gtf_name}.txt
	
			annotatePeaks.pl ${diffbind_dir}/Diffbind_consensus_peak_set_for_Homer.bed \
  			${genome_dir}/${genome}_subset.fa -gtf ${k}.gtf \
  			-CpG > ${homer_dir}/Diffbind_consensus_peak_set_annotated_${gtf_name}.txt
		fi
	done
fi




### generate BigBed file to add peak set to UCSC genome browser session
bedToBigBed ${diffbind_dir}/Diffbind_consensus_peak_set_ucsc_gimmescan.bed \
  ${genome_dir}/${genome}_subset_ucsc.chrom.sizes \
  ${diffbind_dir}/Diffbind_consensus_peak_set_ucsc.bb




### extract gene information from gtf file for downstream R analysis
## No longer necessary
awk 'BEGIN {OFS="\t"} {FS="\t"} {print $9}' ${genome_dir}/${genome}_subset.gtf > genes_1.txt

# subset to genes with ENSDART ID
grep "ENSDART" genes_1.txt > genes_2.txt

# unpack gene names
awk 'BEGIN {OFS="\t"} {FS=";"} {print $1,$3,$6}' genes_2.txt > genes_3.txt

# reduce on gene level
sort genes_3.txt | uniq > genes_4.txt

# get rid of unnecessary columns
awk 'BEGIN {OFS="\t"} {FS="\""} {print $2,$4,$6}' genes_4.txt > genes_5.txt

# get rid of unnecessary rows
grep -v "ensembl_havana" genes_5.txt > ${homer_dir}/genes.txt

rm genes_*.txt




echo All done!
