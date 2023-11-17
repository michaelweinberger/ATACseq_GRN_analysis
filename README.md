# ATAC-seq and gene regulatory network analysis
---
Scripts for an automated analysis pipeline to process bulk ATAC-seq and RNA-seq data, including:
- Mapping of bulk ATAC-seq fastq files via [bowtie](https://bowtie-bio.sourceforge.net/manual.shtml) [1] 
- Mapping of bulk RNA-seq fastq files via [STAR](https://github.com/alexdobin/STAR) [2]
- Peak calling via the Python [MACS3](https://github.com/macs3-project/MACS) [3] package
- Genomic feature annotation of peak regions via [Homer](http://homer.ucsd.edu/homer/ngs/annotation.html) [4]
- Differential peak accessibility analysis using the R [DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) [5] and [DESeq2](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) [6] packages
- Motif calling across peak regions via [gimmescan](https://gimmemotifs.readthedocs.io/en/master/reference.html#command-gimme-scan) [7]
- Gene regulatory network (GRN) construction via the Python [Ananse](https://anansepy.readthedocs.io/en/master/) [8] package



## Usage
---
The pipeline is designed to be run on a high performance cluster, via the Slurm scheduler. Modules loaded in the scripts that we used to run the analysis are:
R-base/4.2.0, python-base/3.10.7, python-base/3.10.7, bowtie/1.2.3, samtools/1.10, bedtools/2.29.2, ucsctools/385, STAR/2.7.10b, subread/2.0.6, homer/20201202, gimmemotifs/20210114

To run the pipeline:
1. Clone the repository via `$ git clone https://github.com/michaelweinberger/ATACseq_GRN_analysis.git`
2. Adjust the `User defined variables` section of the **01_PARENT_script.sh** script:
- `script_dir`   Directory containing all scripts copied from https://github.com/michaelweinberger/ATACseq_GRN_analysis/scripts/
- `out_dir`   Directory containing all output, will be created if non-existent
- `species`   Name of species that sequencing data was generated in, one of "human", "mouse", "chicken" or "zebrafish"
- `peak_width`   The length in bp that each peak region in the consensus peak set will be adjusted to
- `basic_analysis`   Indicates whether mapping, peak calling, Homer annotation, DiffBind analysis and gimmescan motif calling should be performed, one of “Yes” or “No”. 
- `conditions`   List containing the names of the experimental conditions to be compared
- `ATAC_fastq_dirs`   List containing the file paths of the directories containing bulk ATAC-seq fastq files. Each directory should contain all fastq files of a single condition specified in `conditions`. The order of directory entries should match the order of condition listed in `conditions`.
- `RNA_fastq_dirs`   Optional list containing the file paths of the directories containing bulk RNA-seq fastq files. Each directory should contain all fastq files of a single condition specified in `conditions`. The order of directory entries should match the order of condition listed in `conditions`. RNA-seq data will be used to annotate peak regions exclusively to genes expressed in the condition analysed (during DiffBind differential accessibility analysis). If `RNA_fastq_dirs` is empty, peaks will be annotated using all genes in the gtf file.
- `ananse_analysis`   Indicates whether gene regulatory network construction and analysis should be performed, one of “Yes” or “No”. Basic analysis needs to be run before Ananse analysis can be run.
- `conditions_ananse_RNA`   List containing the names of the RNA-seq experimental conditions to be compared
- `conditions_ananse_ATAC`   List containing the names of the ATAC-seq experimental conditions to be compared
- `RNA_fastq_dirs_ananse`   List containing the file paths of the directories containing bulk RNA-seq fastq files to be used during Ananse network construction. The order of directory entries should match the order of condition listed in `conditions_ananse_RNA`. May be the same as `RNA_fastq_dirs`.

3. Finally, start the analysis via `$ sbatch 01_PARENT_script.sh`



## Pipeline details
---
1. **2_1_genome_files_SH.sh**
Generates genome files for the specified species.
- Input: None
- Output: `[out_dir]/genomes` directory containing genome `fasta` and `gtf` files, subset to main chromosomes, as well as `chromosome sizes` and repeat regions `rmsk.txt` file

2. **2_2_loop_mapping_SH.sh**
Generates bowtie [1] genome index, concatenates ATAC-seq fastq files across sequencing lanes, maps concatenated fast files, generates bam files, removes duplicate reads, counts number of entries in final bam files and generates bigwig files (e.g. for UCSC genome browser session).
- Input: `fastq` files located in directories specified in `ATAC_fastq_dirs`, [genome]_subset.fa in `[out_dir]/genomes`
- Output: Bowtie genome index in [out_dir]/genomes; concatenated fastq files in [out_dir]/cat_fastq_ATAC; bowtie.log files in [out_dir]/bowtie_[genome]_mapped; .rmdup.bam, bedgraph (.bg) and bigwig (.bw) files in [out_dir]/bowtie_[genome]_mapped/BAM_files

3. **2_3_loop_peak_calling_SH.sh**
Performs peak calling via MACS [3].
- Input: .rmdup.bam files in [out_dir]/bowtie_[genome]_mapped/BAM_files
- Output: peaks.narrowPeak, control_lambda.bdg, peaks.xls, summits.bed and treat_pileup.bdg files in [out_dir]/MACS

4. **3_DiffBind_input_prep.R**
Removes peaks overlapping repeat regions in out_dir/genomes/rmsk.txt by at least 50% of peak length. Constructs a consensus peak set across all samples using DiffBind [5], adjusts peak length to `peak_width` and counts entries in bam files mapping to consensus peak regions.
- Input: peaks.narrowPeak files in [out_dir]/MACS; .rmdup.bam files in [out_dir]/bowtie_[genome]_mapped/BAM_files; rmsk.txt and .chrom.sizes file in [out_dir]/genomes; conditions.txt file in [out_dir]
- Output: peaks.narrowPeak files with peaks overlapping repetitive elements removed in [out_dir]/MACS/No_repetitive_elements; Diffbind_sample_sheet.csv, Diffbind_consensus_peak_set rds and bed files, sample correlation pdf files in [out_dir]/DiffBind; [genome]_subset_ucsc.chrom.sizes file with chromosome sizes in UCSC nomenclature in [out_dir]/genomes

5. **4_1_loop_mapping_RNA_SH.sh** (optional)
Generates STAR [2] genome index, concatenates RNA-seq fastq files across sequencing lanes, maps concatenated fast files, generates bam files, removes duplicate reads, counts number of entries in final bam files, generates bigwig files (e.g. for UCSC genome browser session) and summarises read counts in bam files via featureCounts [9]. 
- Input: fastq files located in directories specified in `RNA_fastq_dirs`; [genome]_subset.fa, [genome]_subset.gtf, [genome]_subset.chrom.sizes in [out_dir]/genomes
- Output: STAR genome index in [out_dir]/genomes/STAR_[genome]_subset; concatenated fastq files in [out_dir]/cat_fastq_RNA; STAR log files in [out_dir]/STAR_[genome]_mapped/[sample]_mapped; .rmdup.bam, bedgraph (.bg) and bigwig (.bw) files in [out_dir]/STAR_[genome]_mapped/BAM_files; featureCounts output file in [out_dir]/STAR_[genome]_mapped

6. **4_2_featureCounts_output_processing.R** (optional)
Processes raw featureCounts output file.
- Input: featureCounts_raw.txt in [out_dir]/STAR_[genome]_mapped; [genome]_subset.gtf in [out_dir]/genomes
- Output: featureCounts_final.txt in [out_dir]/STAR_[genome]_mapped

7. **4_3_generate_homer_expressed_gene_input.R** (optional)
Generates lists of expressed genes per condition. Considers genes with more than 2.5 FPKM as expressed.
- Input: featureCounts_final.txt in [out_dir]/STAR_[genome]_mapped
- Output: Expressed genes text files in [out_dir]/STAR_[genome]_mapped

8. **5_loop_homer_annotation_SH.sh**
Analyses transcriptional start site (TSS) coverage in bam files, annotates consensus peak set to genomic features using Homer [4]. Optionally subsets gtf file to expressed genes and annotates peaks to expressed genes. Generates bigbed file of consensus peak set (can be added to UCSC genome browser session).
- Input: .rmdup.bam files in [out_dir]/bowtie_[genome]_mapped/BAM_files; [genome]_subset.fa and [genome]_subset.gtf in [out_dir]/genomes; Diffbind_consensus_peak_set_for_Homer.bed and Diffbind_consensus_peak_set_ucsc_gimmescan.bed in [out_dir]/DiffBind; [genome]_subset_ucsc.chrom.sizes file in [out_dir]/genomes; Expressed genes text files in [out_dir]/STAR_[genome]_mapped (optional)
- Output: Tag directory for each sample, Tag_density_tss.txt, Diffbind_consensus_peak_set_annotated text files in [out_dir]/homer_peak_annotation; Diffbind_consensus_peak_set_ucsc.bb bigbed file in [out_dir]/DiffBind; Expressed genes gtf files in [out_dir]/STAR_[genome]_mapped (optional)

9. **6_Homer_annotation.R**
Processes Homer [4] peak annotation output, plots TSS coverage and genomic feature frequencies within consensus peak set.
- Input: Diffbind_consensus_peak_set_annotated text files in [out_dir]/homer_peak_annotation
- Output: Diffbind_annotated_peaks_gene_names rds file, Diffbind_annotated_peaks_gene_names_genomic_feat_freq.pdf, Diffbind_annotated_peaks_gene_names_genomic_feat_gc.pdf, Diffbind_annotated_peaks_gene_names_genomic_feat_cpg.pdf, Diffbind_annotated_peaks_gene_names_tss_enrichment.pdf in [out_dir]/DiffBind

10. **7_DiffBind_DESeq_analysis.R**
Performs differential peak accessibility analysis between conditions using the R DESeq2 [6] package.
- Input: Diffbind_consensus_peak_set, Diffbind_annotated_peaks_gene_names rds files in [out_dir]/DiffBind; conditions.txt in [out_dir]
- Output: DESeq2_ATAC rds files, DESeq2_ATAC csv files with enriched peaks, DESeq2_ATAC volcano plots of logfoldchange versus adjusted p value, DESeq2_ATAC heatmap of accessibility of top differentially accessible peaks in [out_dir]/DiffBind

11. **8_JASPAR2020_gimmescan_SH.sh**
Generates JASPAR2020 [10] transcription factor motif pfm file, calls JASPAR2020 motifs across consensus peaks using gimmescan [7] (false positive rate 0.05).
- Input: Diffbind_consensus_peak_set_ucsc_gimmescan.bed in [out_dir]/DiffBind
- Output: JASPAR2020.pfm in [out_dir], gimme_scan_JASPAR2020_n10000_fpr0_05.txt  in [out_dir]/gimmescan

12. **9_loop_mapping_RNA_SH.sh**
Generates STAR [2] genome index, concatenates RNA-seq fastq files across sequencing lanes, maps concatenated fast files, generates bam files, removes duplicate reads, counts number of entries in final bam files, generates bigwig files (e.g. for UCSC genome browser session) and summarises read counts in bam files via featureCounts [9]. 
- Input: fastq files located in directories specified in `RNA_fastq_dirs_ananse`; [genome]_subset.fa, [genome]_subset.gtf, [genome]_subset.chrom.sizes in [out_dir]/genomes
- Output: STAR genome index in [out_dir]/genomes/STAR_[genome]_subset; concatenated fastq files in [out_dir]/cat_fastq_RNA_ananse; STAR log files in [out_dir]/STAR_[genome]_mapped_ananse/[sample]_mapped; .rmdup.bam, bedgraph (.bg) and bigwig (.bw) files in [out_dir]/STAR_[genome]_mapped_ananse/BAM_files; featureCounts output file in [out_dir]/STAR_[genome]_mapped_ananse

13. **10_featureCounts_output_processing.R **
Processes raw featureCounts output file.
- Input: featureCounts_raw.txt in [out_dir]/STAR_[genome]_mapped_ananse; [genome]_subset.gtf in [out_dir]/genomes
- Output: featureCounts_final.txt in [out_dir]/STAR_[genome]_mapped_ananse

14. **11_Ananse_input_preparation_motif2factors.R**
Generates JASPAR2020 motif2factors file, containing gene names of specified `species`.
- Input: [genome]_subset.gtf in [out_dir]/genomes
- Output: JASPAR2020.motif2factors.txt in [out_dir]

15. **12_Ananse_input_preparation.R**
Generates input for gene regulatory network construction.
- Input: [genome]_subset.gtf in [out_dir]/genomes; featureCounts_final.txt in [out_dir]/STAR_[genome]_mapped_ananse; conditions_ananse_RNA.txt in [out_dir]/Ananse/Ananse_input; 
- Output: ANANSE_TPM.txt files per sample, ANANSE_[genome]_gene_positions.bed in [out_dir]/Ananse/Ananse_input

16. **13_Ananse_SH.sh**
Constructs GRNs using Ananse [8].
- Input: .rmdup.bam files in [out_dir]/STAR_[genome]_mapped_ananse/BAM_files; _peaks.narrowPeak files in [out_dir]/MACS/No_repetitive_elements; [genome]_subset.fa in [out_dir]/genomes; JASPAR2020.pfm in [out_dir]; ANANSE_TPM.txt files per sample, ANANSE_[genome]_gene_positions.bed in [out_dir]/Ananse/Ananse_input
- Output: binding.h5 files in [out_dir]/Ananse/Ananse_output/[condition]_jaspar2020_binding; Network_[condition]_jaspar2020_pro_enh.txt in [out_dir]/Ananse/Ananse_output

17. **14_Ananse_network_follow_up.R**
Analyses number of connections and TFs within GRNs, subsets GRNs (probability above 0.75), generates centrality comparison dot plots between condition, plots centrality heatmap
- Input: JASPAR2020.motif2factors.txt in [out_dir]; conditions_ananse_RNA.txt in [out_dir]/Ananse/Ananse_input; _jaspar2020_pro_enh.txt files in [out_dir]/Ananse/Ananse_output
- Output: Ananse_network subset text files, _connections.pdf,  _factors.pdf, _factor_binding.pdf,  Ananse_network_centrality_jaspar2020_[condition].csv, ananse_network_dot_jaspar2020_ pdf file comparing centrality values, ananse_network_heatmap_jaspar2020.pdf in [out_dir]/Ananse/Ananse_output



## References
---

1.	Langmead, B., Trapnell, C., Pop, M., and Salzberg, S.L. (2009). Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10, R25. 10.1186/gb-2009-10-3-r25.
2.	Dobin, A., Davis, C.A., Schlesinger, F., Drenkow, J., Zaleski, C., Jha, S., Batut, P., Chaisson, M., and Gingeras, T.R. (2013). STAR: ultrafast universal RNA-seq aligner. Bioinformatics 29, 15-21. 10.1093/bioinformatics/bts635.
3.	Zhang, Y., Liu, T., Meyer, C.A., Eeckhoute, J., Johnson, D.S., Bernstein, B.E., Nusbaum, C., Myers, R.M., Brown, M., Li, W., and Liu, X.S. (2008). Model-based analysis of ChIP-Seq (MACS). Genome Biol 9, R137. 10.1186/gb-2008-9-9-r137.
4.	Heinz, S., Benner, C., Spann, N., Bertolino, E., Lin, Y.C., Laslo, P., Cheng, J.X., Murre, C., Singh, H., and Glass, C.K. (2010). Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. Mol Cell 38, 576-589. 10.1016/j.molcel.2010.05.004.
5.	Stark, R., and Brown, G. (2011). DiffBind: differential binding analysis of ChIP-Seq peak data. .
6.	Love, M.I., Huber, W., and Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15, 550. 10.1186/s13059-014-0550-8.
7.	van Heeringen, S.J., and Veenstra, G.J. (2011). GimmeMotifs: a de novo motif prediction pipeline for ChIP-sequencing experiments. Bioinformatics 27, 270-271. 10.1093/bioinformatics/btq636.
8.	Xu, Q., Georgiou, G., Frolich, S., van der Sande, M., Veenstra, G.J.C., Zhou, H., and van Heeringen, S.J. (2021). ANANSE: an enhancer network-based computational approach for predicting key transcription factors in cell fate determination. Nucleic Acids Res. 10.1093/nar/gkab598.
9.	Liao, Y., Smyth, G.K., and Shi, W. (2014). featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30, 923-930. 10.1093/bioinformatics/btt656.
10.	Fornes, O., Castro-Mondragon, J.A., Khan, A., van der Lee, R., Zhang, X., Richmond, P.A., Modi, B.P., Correard, S., Gheorghe, M., Baranasic, D., et al. (2020). JASPAR 2020: update of the open-access database of transcription factor binding profiles. Nucleic Acids Res 48, D87-D92. 10.1093/nar/gkz1001.


