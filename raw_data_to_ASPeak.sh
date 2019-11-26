# !/bin/bash

### Maria Mikedis mikedis@wi.mit.edu
### 4/12/18

DIR="/lab/solexa_page/maria/iCLIP_DAZL_lepts"
cd ${DIR}



##############
### check quality report on raw data
##############

cd 170921_iCLIP_rawdata/170914_WIGTC-HISEQ2B_CBFT3ANXX/QualityScore/
bsub -K "fastqc unknown-s_3_1_sequence.txt.gz"
wait



##############
### quality trim reads
##############
###	require minimum quality score of 20, minimum length of 24

cd ${DIR}
mkdir 170921_quality_trimmed
cd 170921_quality_trimmed
ln -s /lab/solexa_page/maria/iCLIP_DAZL_lepts/170921_iCLIP_rawdata/170914_WIGTC-HISEQ2B_CBFT3ANXX/QualityScore/unknown-s_3_1_sequence.txt.gz .

bsub -K "cutadapt -q 20 -m 24 -o unknown-s_3_1_sequence.trimmed.txt.gz unknown-s_3_1_sequence.txt.gz"
wait


##############
### remove PCR duplicates
##############
###	require minimum quality score of 20, minimum length of 24

cd ${DIR}
mkdir 170921_removal_PCR_duplicates
cd 170921_removal_PCR_duplicates

cp /lab/solexa_page/maria/iCLIP_DAZL_lepts/170921_quality_trimmed/unknown-s_3_1_sequence.trimmed.txt.gz . 
bsub -K gzip -d unknown-s_3_1_sequence.trimmed.txt.gz
wait

bsub -K "fastx_collapser -v -i  unknown-s_3_1_sequence.trimmed.txt > unknown-s_3_1_sequence.trimmed_PCR.dup.removed.txt"
wait 
rm -f unknown-s_3_1_sequence.trimmed.txt




##############
### demultiplex
##############

cd ${DIR}
mkdir 170921_demultiplex
cd 170921_demultiplex

INPUT="/lab/solexa_page/maria/iCLIP_DAZL_lepts/170921_removal_PCR_duplicates/unknown-s_3_1_sequence.trimmed_PCR.dup.removed.txt"

###remove first three nucleotides (random barcode)
fastx_trimmer -f 4 -i ${INPUT} -o unknown-s_3_1_sequence.trimmed_PCR.dup.removed_ran.bc.trim.txt

cat unknown-s_3_1_sequence.trimmed_PCR.dup.removed_ran.bc.trim.txt | fastx_barcode_splitter.pl --bcfile barcodes --prefix demultiplexed_ --bol 

bsub "gzip unknown-s_3_1_sequence.trimmed_PCR.dup.removed_ran.bc.trim.txt"

### trim off regular barcodes (4 nucleotides) and 2 nucleotides of random barcode on  5' end
for FILE in `ls demultiplexed* | sed s/@//`; do
	fastx_trimmer -f 7 -i ${FILE} -o ${FILE}_trimmed
done


##############
### map reads and identify crosslinked nucleotides
##############

cd ${DIR}
now1="$(date +'%Y%m%d')"
mkdir ${now1}_iCLIP_STAR
CLIPDIR="/lab/solexa_page/maria/iCLIP_DAZL_lepts/${now1}_iCLIP_STAR"
cd ${CLIPDIR}

ln -s /lab/solexa_page/maria/iCLIP_DAZL_lepts/170921_demultiplex/*trimmed .


### run STAR for iCLIP reads wtih 30nt overhang
### modify output bam files so only crosslinked nt is represented by each read
### index bam files

INDEX="/lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/STAR/STAR_output_overhang30/"
jobnum=1

for FILE in demultiplexed_MM247_1_DAZL_CLIP_trimmed demultiplexed_MM247_2_IgG_CLIP_trimmed demultiplexed_MM253_1_1_DAZL_CLIP_trimmed demultiplexed_MM253_1_2_IgG_CLIP_trimmed demultiplexed_MM253_2_1_DAZL_CLIP_trimmed demultiplexed_MM253_2_2_IgG_CLIP_trimmed; do 

	bsub -K -J "iCLIP.jobs${jobnum}" -n8 -R "span[hosts=1]" "/nfs/apps/STAR/2.5.4b/STAR  --genomeDir ${INDEX} --readFilesIn ${FILE}  --outFileNamePrefix ${FILE}_STAR_ --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --outFilterMismatchNmax 2 --runThreadN 8 --outSAMattributes None --outReadsUnmapped Fastx"
	wait
	
	bsub -K  -J "iCLIP.jobs${jobnum}.1" "python converting_mapped_reads_to_crosslinked_sites.py ${FILE}_STAR_Aligned.out.sam ${FILE}_STAR_Aligned.out_CL.site.only.sam"
	wait

	bsub -J "iCLIP.jobs${jobnum}.2" "samtools view -bT /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10.fa  ${FILE}_STAR_Aligned.out.sam >|  ${FILE}_STAR_Aligned.out.bam"
	bsub -K  -J "iCLIP.jobs${jobnum}.3" "samtools view -bT /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10.fa  ${FILE}_STAR_Aligned.out_CL.site.only.sam >|  ${FILE}_STAR_Aligned.out_CL.site.only.bam"
	wait
	bsub  -J "iCLIP.jobs${jobnum}.4" "rm -f ${FILE}_STAR_Aligned.out.sam"
	bsub -J "iCLIP.jobs${jobnum}.5" "rm -f ${FILE}_STAR_Aligned.out_CL.site.only.sam"

	bsub  -K -J "iCLIP.jobs${jobnum}.6" "samtools sort ${FILE}_STAR_Aligned.out.bam >| ${FILE}_STAR_Aligned.out_sorted.bam"
	wait
	bsub  -K  -J "iCLIP.jobs${jobnum}.7" "samtools sort ${FILE}_STAR_Aligned.out_CL.site.only.bam >| ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam" 
	wait

	bsub  -J "iCLIP.jobs${jobnum}.8" "samtools index ${FILE}_STAR_Aligned.out_sorted.bam ${FILE}_STAR_Aligned.out_sorted.bam.bai"
	bsub  -J "iCLIP.jobs${jobnum}.9" "samtools index ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam.bai"

	jobnum=$((jobnum + 1))

done



############
### Call peaks
############
### Run ASPeak, using crosslinked nucletotides only
###				using IgG as background
### 		    using ASPeak bed files where overlapping genomic regions have already been removed based on this hierarchy: 
###				3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### Analyze library prep batches separately (MM247 samples represent library batch 1; MM253 samples represent library batch 2)



cd ${DIR}
now3="$(date +'%Y%m%d')"
mkdir ${now3}_ASPeak
cd ${now3}_ASPeak


### link files containing crosslinked sites
for FILE in demultiplexed_MM247_1_DAZL_CLIP_trimmed demultiplexed_MM247_2_IgG_CLIP_trimmed demultiplexed_MM247_3_DAZL_noCL_trimmed demultiplexed_MM247_4_input_trimmed demultiplexed_MM253_1_1_DAZL_CLIP_trimmed demultiplexed_MM253_1_2_IgG_CLIP_trimmed demultiplexed_MM253_1_3_DAZL_noCL_trimmed demultiplexed_MM253_1_4_input_trimmed demultiplexed_MM253_2_1_DAZL_CLIP_trimmed demultiplexed_MM253_2_2_IgG_CLIP_trimmed demultiplexed_MM253_2_3_DAZL_noCL_trimmed demultiplexed_MM253_2_4_input_trimmed; do 
bsub "ln -s ${CLIPDIR}/${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam ."
done

CLIP1=demultiplexed_MM247_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
CLIP2=demultiplexed_MM253_1_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
CLIP3=demultiplexed_MM253_2_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG1=demultiplexed_MM247_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG2=demultiplexed_MM253_1_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG3=demultiplexed_MM253_2_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam


### Run ASPeak
bsub "ASPeak.pl -param ~/bin/ASPeak/scripts/default_parameters_mm10.txt -gapnumber 0 -bed /lab/solexa_page/maria/genome_files/ASPeak_171122/bed_files_modified_for_ASPeak_overlapping_regions_excluded -lib ${CLIP2}:${CLIP3} -outdir ASPeak_IgG_controls_gapno0_nornaseq_MM253_only_hierarchy_applied -control ${IgG2}:${IgG3}"
bsub "ASPeak.pl -param ~/bin/ASPeak/scripts/default_parameters_mm10.txt -gapnumber 0 -bed /lab/solexa_page/maria/genome_files/ASPeak_171122/bed_files_modified_for_ASPeak_overlapping_regions_excluded -lib ${CLIP1} -outdir ASPeak_IgG_controls_gapno0_nornaseq_MM247_only_hierarchy_applied -control ${IgG1}"

