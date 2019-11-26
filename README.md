# Gura_et_al_2019_PLOS_Genetics
Custom scripts used to analyze DAZL iCLIP data from Gura et al., 2019 PLOS Genetics

raw_data_to_ASPeak.sh
  Calls iCLIP peaks from raw sequencing data

converting_mapped_reads_to_crosslinked_sites
  Called by raw_data_to_ASpeak.sh to convert mapped reads to crosslinked sites

post_ASPeak_analysis.sh
Analyzes peaks called by ASPeak

extract_3UTR_sequences.R
  Called by post_ASpeak_analysis.sh
  Extracts 3' UTR sequences from mm10 genome using UCSC table "refGene"

make_venn_diagram_with_lines.R
  Called by post_ASpeak_analysis.sh
  Creates a venn diagram to show overlap of peaks called among 3 biological replicates
  
obtain_3UTRseq_from_bed_190220.R
  Called by post_ASpeak_analysis.sh
  Takes a bed file, list of Refseq ids, and number of nucleotides by which to expand the bed file sequence along the transcript
  Provides 3' UTR sequence
