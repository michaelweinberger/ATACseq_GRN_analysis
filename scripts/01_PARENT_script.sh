#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=2-00:00:00
#SBATCH --mem=60G
#SBATCH --cpus-per-task=12




#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

#################################################################################################################################################


# specify directory containing parent and child scripts
script_dir="/ceph/project/tsslab/mweinber/tmp/scripts"


# specify output directory
# if directory does not exist, it will be generated
out_dir="/ceph/project/tsslab/mweinber/tmp/out"


# specify species ("human", "mouse", "chicken" or "zebrafish")
species="zebrafish"


# specify desired length of peak regions in consensus peak set [bp]
peak_width=500



##############   Mapping + DiffBind analysis   ##############

# indicate if basic ATAC-seq analysis should be run ("Yes" or "No")
basic_analysis="Yes"


# specify the names of the different experimental conditions
conditions=("tcf21_larval"
 "tcf21_cryo"
)


# specify the location of directories containing bulk ATAC-seq fastq files
# all fastq files belonging to samples of the same condition should be in the SAME directory, each condition should have a SEPARATE directory
# -> all fastq files in a specified directory will be analysed as belonging to the respective condition (specified above)
ATAC_fastq_dirs=("/ceph/project/tsslab/mweinber/tmp/ATAC_tcf21_larval"
 "/ceph/project/tsslab/mweinber/tmp/ATAC_tcf21_cryo"
)


# OPTIONAL: specify the location of directories containing bulk RNA-seq fastq files
# the order of directory entries in "RNA_fastq_dirs" should match the order of corresponding directories in "ATAC_fastq_dirs"
# this is used to annotate peak regions only to genes expressed in the individual conditions, instead of annotating to all genes in the genome gtf file
# for each condition, the respective expressed-genes-only annotation will then be used in the downstream DiffBind analysis
# if "RNA_fastq_dirs" is empty, peak regions will be annotated to all genes in the genome gtf file
RNA_fastq_dirs=("/ceph/project/tsslab/mweinber/tmp/RNA_tcf21_larval"
 "/ceph/project/tsslab/mweinber/tmp/RNA_tcf21_cryo"
)



##############   ANANSE analysis   ##############

# indicate if gene regulatory network (GRN) construction + analysis should be run ("Yes" or "No")
ananse_analysis="Yes"


# specify directories containing bulk RNA-seq fastq files (needed for Ananse GRN construction)
# all fastq files belonging to samples of the same condition should be in the SAME directory, each condition should have a SEPARATE directory
# fastq file names in each directory should contain the corresponding condition name (see below)
RNA_fastq_dirs_ananse=("/ceph/project/tsslab/mweinber/tmp/RNA_tcf21_larval"
 "/ceph/project/tsslab/mweinber/tmp/RNA_tcf21_cryo"
)


# specify conditions to be analysed with Ananse
# the order of entries in "conditions_ananse_RNA" should match the order of directory entries in "RNA_fastq_dirs_ananse"
# conditions with matching positions in "conditions_ananse_RNA" and "conditions_ananse_ATAC" will be jointly used for GRN construction
conditions_ananse_RNA=("tcf21_larval"
 "tcf21_cryo"
)

conditions_ananse_ATAC=("tcf21_larval"
 "tcf21_cryo"
)


#################################################################################################################################################




### Analysis

# load modules
module load R-cbrg


# specify additional genome information
if test "$species" = "human" ; then
	species_latin="homo_sapiens"
	species_latin_2="Homo_sapiens"
	genome="GRCh38"
	genome_ucsc="hg38"
	chrom_number=22 #X,Y chromosomes excluded
elif test "$species" = "mouse" ; then
	species_latin="mus_musculus"
	species_latin_2="Mus_musculus"
	genome="GRCm39"
	genome_ucsc="mm39"
	chrom_number=19 #X,Y chromosomes excluded
elif test "$species" = "chicken" ; then
	species_latin="gallus_gallus"
	species_latin_2="Gallus_gallus"
	genome="bGalGal1.mat.broiler.GRCg7b"
	genome_ucsc="galGal6"
	chrom_number=33 #W,Z chromosomes excluded
elif test "$species" = "zebrafish" ; then
	species_latin="danio_rerio"
	species_latin_2="Danio_rerio"
	genome="GRCz11"
	genome_ucsc="danRer11"
	chrom_number=25
fi


# create output directory if it does not exist
[ ! -d "$out_dir" ] && mkdir -p "$out_dir"


# specify location of directory containing genome files
genome_dir=${out_dir}/genomes



## run mapping, peak calling, DiffBind, gimmescan..
if test $basic_analysis = "Yes"; then

	# save conditions array, to be used in R scripts
	# replace dashes and dots in condition names with underscores
	num_cond=${#conditions[@]}
	for (( i=0; i<$num_cond; i++ ))
	do
		if test $i = 0;then
			echo ${conditions[i]} | tr '-' '_' | tr '.' '_' > ${out_dir}/conditions.txt
		else
			echo ${conditions[i]} | tr '-' '_' | tr '.' '_' >> ${out_dir}/conditions.txt
		fi
	done

	# run genome file preparation script
	source ${script_dir}/2_1_genome_files_SH.sh

	# run mapping script
	source ${script_dir}/2_2_loop_mapping_SH.sh

	# run peak calling script
	source ${script_dir}/2_3_loop_peak_calling_SH.sh

	# run DiffBind consensus peak construction 
	# create directory for DiffBind output
	diffbind_dir="${out_dir}/DiffBind"
	[ ! -d "$diffbind_dir" ] && mkdir -p "$diffbind_dir"

	R -f ${script_dir}/3_DiffBind_input_prep.R --args out_dir=${diffbind_dir} peak_dir=${peak_dir} bam_dir=${bam_dir} \
  	conditions_dir=${out_dir} genome_dir=${genome_dir} genome=${genome} chr_sizes="${genome_dir}/${genome}.chrom.sizes" \
  	repeats_file="${genome_dir}/rmsk.txt" chrom_number=${chrom_number} peak_width=${peak_width}


	if [ ${#RNA_fastq_dirs[@]} > 0 ]; then
		# run bulk RNAseq mapping script
		source ${script_dir}/4_1_loop_mapping_RNA_SH.sh

		# run featureCounts output modification script
		if [ ! -f "${sam_dir_2}/featureCounts_final.txt" ]; then
			R -f ${script_dir}/4_2_featureCounts_output_processing.R --args out_dir=${sam_dir_2} \
			feature_counts_file_raw=${sam_dir_2}/featureCounts_raw.txt gtf_file="${genome_dir}/${genome}_subset.gtf"
		fi

                # create list of expressed genes for each condition 
		R -f ${script_dir}/4_3_generate_homer_expressed_gene_input.R --args out_dir=${sam_dir_2} \
		conditions_dir=${out_dir} feature_counts_file=${sam_dir_2}/featureCounts_final.txt
	fi

	# run Homer annotation script
	source ${script_dir}/5_loop_homer_annotation_SH.sh

	# run Homer downstream analysis script
	R -f ${script_dir}/6_Homer_annotation.R --args out_dir=${diffbind_dir} conditions_dir=${out_dir} homer_dir=${homer_dir} \
	genome=${genome}

	# run DESeq2 differential peak accessibility analysis script
	R -f ${script_dir}/7_DiffBind_DESeq_analysis.R --args out_dir=${diffbind_dir} conditions_dir=${out_dir} 

	# avoid conflict with gimmescan Python requirement
	module unload python-cbrg

	# run gimmescan motif calling script
	source ${script_dir}/8_JASPAR2020_gimmescan_SH.sh

	# avoid conflict with Ananse
	module unload python-base/3.8.3
fi



## run Ananse analysis
if test $ananse_analysis = "Yes"; then

	# define important directories from previous analysis part
	sam_dir="${out_dir}/bowtie_${genome}_mapped"
	bam_dir="${sam_dir}/BAM_files"
	peak_dir="${out_dir}/MACS"

	# create Ananse directories
	ananse_input_dir="${out_dir}/Ananse/Ananse_input"
	[ ! -d "$ananse_input_dir" ] && mkdir -p "$ananse_input_dir"
	ananse_output_dir="${out_dir}/Ananse/Ananse_output"
	[ ! -d "$ananse_output_dir" ] && mkdir -p "$ananse_output_dir"

	# save conditions array, to be used in R scripts
	# replace dashes and dots in condition names with underscores
	num_cond=${#conditions_ananse_ATAC[@]}
	for (( i=0; i<$num_cond; i++ ))
	do
		if test $i = 0;then
			echo ${conditions_ananse_ATAC[i]} | tr '-' '_' | tr '.' '_' > ${ananse_input_dir}/conditions_ananse_ATAC.txt
		else
			echo ${conditions_ananse_ATAC[i]} | tr '-' '_' | tr '.' '_' >> ${ananse_input_dir}/conditions_ananse_ATAC.txt
		fi
	done

	num_cond=${#conditions_ananse_RNA[@]}
	for (( i=0; i<$num_cond; i++ ))
	do
		if test $i = 0;then
			echo ${conditions_ananse_RNA[i]} | tr '-' '_' | tr '.' '_' > ${ananse_input_dir}/conditions_ananse_RNA.txt
		else
			echo ${conditions_ananse_RNA[i]} | tr '-' '_' | tr '.' '_' >> ${ananse_input_dir}/conditions_ananse_RNA.txt
		fi
	done

	# run bulk RNAseq mapping script
	source ${script_dir}/9_loop_mapping_RNA_SH.sh

	# run featureCounts output modification script
	if [ ! -f "${sam_dir_2}/featureCounts_final.txt" ]; then
		R -f ${script_dir}/10_featureCounts_output_processing.R --args out_dir=${sam_dir_2} \
		feature_counts_file_raw=${sam_dir_2}/featureCounts_raw.txt gtf_file="${genome_dir}/${genome}_subset.gtf"
	fi

	# run Ananse input preparation Jaspar2020 motif2factors file generation script
	if [ ! -f "${out_dir}/JASPAR2020.motif2factors.txt" ]; then
		R -f ${script_dir}/11_Ananse_input_preparation_motif2factors.R --args out_dir=${out_dir} gtf="${genome_dir}/${genome}_subset.gtf" \
		species=${species}
	fi

	# run Ananse input preparation script
	R -f ${script_dir}/12_Ananse_input_preparation.R --args peak_dir=${peak_dir} bam_dir=${bam_dir} out_dir=${ananse_input_dir} \
        peak_width=${peak_width} chr_sizes="${genome_dir}/${genome}_subset.chrom.sizes" chrom_number=${chrom_number} \
        genome=${genome} species=${species} gtf_file="${genome_dir}/${genome}_subset.gtf" \
	feature_counts_file=${sam_dir_2}/featureCounts_final.txt

	# run Ananse
	source ${script_dir}/13_Ananse_SH.sh

	# run Ananse output processing script
	R -f ${script_dir}/14_Ananse_network_follow_up.R --args ananse_input_dir=${ananse_input_dir} ananse_output_dir=${ananse_output_dir} \
        peak_width=${peak_width} chrom_number=${chrom_number} out_dir=${out_dir}
fi


echo All done!
