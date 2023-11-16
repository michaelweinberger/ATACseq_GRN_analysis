#!/bin/sh

#Format of --time is DAYS-HOURS:MINUTES:SECONDS
#SBATCH --time=0-06:00:00
#SBATCH --mem=10G



#############################################################
####                                                    #####
####              User defined variables                #####
####                                                    #####
#############################################################

########################################################################################################################


########################################################################################################################



module load gimmemotifs/20210114
# This command also loads requirement: python-base/3.8.3 bedtools/2.29.2 ucsctools/385 R-base/4.2.0 gsl/2.6 hdf5/1.10.7 gdal/3.4.1 sqlite/3.37.2 proj/8.1.1 geos/3.9.0 jags/4.3.0 java/17.0.1 cmake/3.17.2 R-cbrg/current homer/20201202 perl/5.32.0 meme/5.4.1



### create directory
gimmescan_dir="${out_dir}/gimmescan"
[ ! -d "$gimmescan_dir" ] && mkdir -p "$gimmescan_dir"



cd ${out_dir}



### process JASPAR2020 pfm transcription factor motif file, downloaded from JASPAR2020 Core Vertebrates in MEME format as .txt file
echo Generating JASPAR2020 motif file
curl https://jaspar.elixir.no/download/data/2020/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt

# extract lines that start with MOTIF into new file
grep '^MOTIF' JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt > JASPAR2020_motifs_extracted.txt

# remove first 9 lines of the files, remove lines beginning with URL + following line, remove lines beginning with "letter"
cat JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt | sed -e '1,9d' | sed '/^URL/,+1 d' | sed '/^letter/ d' > JASPAR2020_raw.pfm 

# replace MOTIF at begin of line with >
sed -i 's/MOTIF />/g' JASPAR2020_raw.pfm

# replace double whitespaces with tab
sed -i 's/  /\t/g' JASPAR2020_raw.pfm

# delete leading whitespace
sed -i 's/^ //g' JASPAR2020_raw.pfm

# delete substrings following whitespace
awk -F'[ ]' '{print $1}' < JASPAR2020_raw.pfm > ${out_dir}/JASPAR2020.pfm

rm JASPAR2020_raw.pfm
rm JASPAR2020_CORE_vertebrates_non-redundant_pfms_meme.txt



### call transcription factor binding sites across all peaks in the DiffBind peak set
# install genome for use with gimmescan
echo Installing ${genome_ucsc} for gimmemotifs
genomepy install ${genome_ucsc} --provider UCSC --annotation 

wait

# run gimmescan command
if [ ! -f "${gimmescan_dir}/gimme_scan_JASPAR2020_n10000_fpr0_05.txt" ]; then
	echo Running gimmescan analysis using ${genome_ucsc} genome and ${out_dir}/JASPAR2020.pfm on ${diffbind_dir}/Diffbind_consensus_peak_set_ucsc_gimmescan.bed
	gimme scan ${diffbind_dir}/Diffbind_consensus_peak_set_ucsc_gimmescan.bed \
	-g ${genome_ucsc} \
	-p ${out_dir}/JASPAR2020.pfm \
	-n 10000 -f 0.05 -t > ${gimmescan_dir}/gimme_scan_JASPAR2020_n10000_fpr0_05.txt
fi



echo All done!

