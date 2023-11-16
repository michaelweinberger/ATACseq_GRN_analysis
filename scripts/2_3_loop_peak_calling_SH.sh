#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-12:00:00
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



### perform peak calling
peak_dir="${out_dir}/MACS"
[ ! -d "$peak_dir" ] && mkdir -p "$peak_dir"

cd $bam_dir

for j in $(ls *.rmdup.bam | rev | cut -c 11- | rev | uniq)
do
	if [ ! -f "${peak_dir}/${j}_peaks.narrowPeak" ]; then
		# sort BAM files
		samtools sort -n -o ${j}.rmdup.sorted.bam ${j}.rmdup.bam

		# perform peak calling with MACSv3:
		echo Calling peaks in ${j} BAM file
		macs3 callpeak -B -t ${j}.rmdup.sorted.bam -f BAMPE --outdir ${peak_dir} -n ${j} -g ${genome_size} --verbose 1 --call-summits

		# remove unnecessary files
		rm ${j}.rmdup.sorted.bam
	fi
done

wait



echo All done! 

