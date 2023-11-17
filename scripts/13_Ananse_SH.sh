#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=12



# to install ananse
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda config --set use_only_tar_bz2 True
#conda create -n ananse ananse
#conda activate ananse

module load bedtools



### run Ananse analysis, using JASPAR 2020 motif database
# loop through list of conditions
len=${#conditions_ananse_ATAC[@]}

for (( k=0; k<${len}; k++ ))
do
	# TF binding analysis
	if [ ! -f "${ananse_output_dir}/${conditions_ananse_ATAC[k]}_jaspar2020_binding/binding.h5" ]; then

		echo Ananse Binding analysis of ${conditions_ananse_ATAC[k]} using ${ananse_input_dir}/JASPAR2020.pfm database

		bam_input=$(find ${bam_dir} -maxdepth 1 -regex ${bam_dir}/${conditions_ananse_ATAC[k]}.*rmdup\.bam$)
		peak_input=$(find ${peak_dir} -maxdepth 1 -regex ${peak_dir}/.*${conditions_ananse_ATAC[k]}.*_peaks\.narrowPeak)
		echo Bam files: ${bam_input}
		echo Peak files: ${peak_input}
	
		ananse binding -n 12 \
		-A ${bam_input} \
		-r ${peak_input} \
		-o ${ananse_output_dir}/${conditions_ananse_ATAC[k]}_jaspar2020_binding \
		-g ${genome_dir}/${genome}_subset.fa \
		-p ${out_dir}/JASPAR2020.pfm
	fi
	
	wait
	
	# network construction
	if [ ! -f "${ananse_output_dir}/Network_${conditions_ananse_RNA[k]}_jaspar2020_pro_enh.txt" ]; then

		echo Ananse Network construction of ${conditions_ananse_ATAC[k]} and ${conditions_ananse_RNA[k]}

		tpm_input=$(find ${ananse_input_dir} -maxdepth 1 -regex ${ananse_input_dir}/.*${conditions_ananse_RNA[k]}.*_ANANSE_TPM\.txt)
		echo TPM input files: ${tpm_input}

		ananse network -n 12 \
		-e ${tpm_input} \
		-o ${ananse_output_dir}/Network_${conditions_ananse_RNA[k]}_jaspar2020_pro_enh.txt \
		-g ${genome_dir}/${genome}_subset.fa \
		-a ${ananse_input_dir}/ANANSE_${genome}_gene_positions.bed \
		${ananse_output_dir}/${conditions_ananse_ATAC[k]}_jaspar2020_binding/binding.h5 \
		--include-promoter --include-enhancer
	fi
done

wait



echo All done!
