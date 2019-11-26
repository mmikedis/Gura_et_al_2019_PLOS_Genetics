# !/bin/bash

### Maria Mikedis mikedis@wi.mit.edu

### for each analysis, update variable DIR, SAMPLE1, SAMPLE2, SAMPLE3, TPMS, EXPRESSED, FA, ANNO, BAMDIR
### This script will create directory <${ASPEAK_output_directory}/foundpeaks/merged.peaks.191002> and put files into region-specific directories


DIR1="/lab/solexa_page/maria/iCLIP_DAZL_lepts/20191001_ASPeak/ASPeak_IgG_controls_gapno0_nornaseq_MM247_only_hierarchy_applied"
DIR2="/lab/solexa_page/maria/iCLIP_DAZL_lepts/20191001_ASPeak/ASPeak_IgG_controls_gapno0_nornaseq_MM253_only_hierarchy_applied"


RDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/R_scripts_post_ASPeak_analysis"

FA="/nfs/genomes/mouse_mm10_dec_11_no_random/fasta_whole_genome/mm10.fa"
	### fasta of chromosome sequences (each chromosome is a fasta entry)

ANNO="/lab/solexa_page/maria/genome_files/kallisto/for_ASPeak_analysis/mm10_GRCm38_refGene_mRNA_retrogenes_names.txt"
	### annotation file of transcript IDs (field 1) and gene IDs (field 2)

BAMDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180416_iCLIP_STAR"
	### directory containing bam files of aligned CLIP data

SAMPLE1="MM247_1"
SAMPLE2="MM253_1_1"
SAMPLE3="MM253_2_1"


IGG1="MM247_2"
IGG2="MM253_1_2"
IGG3="MM253_2_2"


cd ${DIR1}/foundpeaks
mkdir merged.peaks.191002


### convert peak files to bed files; merge peaks so that one peak can represent multiple transcripts; link to merged.peak directory
### filter peaks for p val <0.0001
### if no peaks meet the designated cut offs, bedtools merge will give the error message: "ERROR: Requested column 4, but database file stdin only has fields 1 - 0."
### bedtools output the strand in the wrong column so modifying bedtools output to follow bed format
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only 
	do
	mkdir ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/
	PEAKS="${DIR1}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed .

done

### for peaks that are found in multiple "regions", remove peaks based on the following hierarchy: 3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### for example, if a peak shows up in both 3' UTR and ncRNA, the hierarchical priortization removes it from ncRNA and leaves it in 3' UTR
### if there is only partial overlap in the peaks, the full peak will be removed

for REGION in intron_coding.only; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE1}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in retrogenes; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done

for REGION in ncRNA.no.rRNA; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done

for REGION in cds; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done



for REGION in utr5_coding.only; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done

	### saving utr3_coding.only peak file with hierachry.applied in name to keep file names consistent; no utr3 peaks are removed
for REGION in utr3_coding.only tRNAs rRNA.only intergenic; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	cat ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done





### repeat for MM253 samples

cd ${DIR2}/foundpeaks
mkdir merged.peaks.191002

### convert peak files to bed files; merge peaks so that one peak can represent multiple transcripts; link to merged.peak directory
### filter peaks for p val <0.0001
### if no peaks meet the designated cut offs, bedtools merge will give the error message: "ERROR: Requested column 4, but database file stdin only has fields 1 - 0."
### bedtools output the strand in the wrong column so modifying bedtools output to follow bed format
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only 
	do
	mkdir "${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/"
	PEAKS="${DIR2}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed .
	PEAKS="${DIR2}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed .
done

### for peaks that are found in multiple "regions", remove peaks based on the following hierarchy: 3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### for example, if a peak shows up in both 3' UTR and ncRNA, the hierarchical priortization removes it from ncRNA and leaves it in 3' UTR
### if there is only partial overlap in the peaks, the full peak will be removed

for REGION in intron_coding.only; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE2}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE3}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in retrogenes; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >|  ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in ncRNA.no.rRNA; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in cds; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}


	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >|  ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done



for REGION in utr5_coding.only; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done

	### saving utr3_coding.only peak file with hierachry.applied in name to keep file names consistent; no utr3 peaks are removed
for REGION in utr3_coding.only tRNAs rRNA.only intergenic; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}


	cat ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	cat ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done





cd ${DIR1}
cd ..
mkdir ASPeak_MM247_MM253_peaks_called_independently_hierarchy_applied
MERGED=/lab/solexa_page/maria/iCLIP_DAZL_lepts/20191001_ASPeak/ASPeak_MM247_MM253_peaks_called_independently_hierarchy_applied
cd ${MERGED}

### identify peaks that are present in at least 2 replicates using bedtools; one loop for regions where hierachy was applied; second loop for other regions
for REGION in cds intron_coding.only ncRNA.no.rRNA retrogenes utr3_coding.only utr5_coding.only; do
	cd ${MERGED}
	mkdir ${REGION}
	cd ${REGION}
	CLIP1=${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP2=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP3=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	ln -s ${CLIP1};
	ln -s ${CLIP2};
	ln -s ${CLIP3};
	bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
done

for REGION in intergenic rRNA.only tRNAs; do
	cd ${MERGED}
	mkdir ${REGION}
	cd ${REGION}
	CLIP1=${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP2=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP3=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	ln -s ${CLIP1};
	ln -s ${CLIP2};
	ln -s ${CLIP3};
	bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
done


### create a master bed file with genomic position information only for peaks that are present in at least 2 of the 3 biological replicates
### I used the following strategy instead of bedtools because bedtools merge was doing funny things when I tried to merge intersected files
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
cd ${MERGED}/${REGION}/
cut -f1-3,6 ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq >| tmp1.tmp
cut -f1-3,6 ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed| sort -k1,1 -k2,2n | uniq >| tmp2.tmp
cut -f1-3,6 ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq  >| tmp3.tmp
cat tmp1.tmp tmp2.tmp >| tmp.tmp
cat tmp.tmp tmp3.tmp | sort -k1,1 -k2,2n | uniq | awk '{FS="\t"; print $1"\t"$2"\t"$3"\t""name""\t""score""\t"$4}'>|  ${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed
rm -f tmp*
done


####pull out replicated peaks from each replicates peak.bed file 
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	PEAKS=${MERGED}/${REGION}/${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed;
	REP1="${DIR1}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP1};
	bedtools intersect -header -s -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	REP2="${DIR2}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP2};
	bedtools intersect -header -s -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	REP3="${DIR2}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP3};
	bedtools intersect -header -s -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
done

### combine bare-bones master bed file with peak information to create a master bed file with more complete info

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}	
	FILE1="${DIR1}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE2="${DIR2}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE3="${DIR3}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	cat ${FILE1} ${FILE2} >| tmp1.tmp
	cat tmp1.tmp ${FILE3} | sort -k1,1 -k2,2n | bedtools merge -header -s -c 4,5,6,7,8,9,10,11,12,13 -o distinct,max,distinct,max,max,min,max,max,max,max -i stdin | uniq |awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' >| ${REGION}.peaks_present.in.at.least.2.replicates.bed; # uniq command gets rid of repeated headers; bedtools behaving weird and output 2 strand columns
	rm -f tmp*;
done



### create a list of genes with peaks, without header

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


for REGION in retrogenes; do 
	cd ${MERGED}/${REGION}
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_exon"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


### calculate overlap in peaks among replicates; produce venn diagrams
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	area1="$(grep -v "#" 	${MERGED}/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq | wc | awk '{print $1}')"
	area2="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	area3="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	n12="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n23="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n13="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n123="$(grep -v "#" ${MERGED}/${REGION}/${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	echo "Replicate No_peaks" >| tmp.tmp
	echo "${SAMPLE1} ${area1}" >> tmp.tmp
	echo "${SAMPLE2} ${area2}">> tmp.tmp
	echo "${SAMPLE3} ${area3}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2} ${n12}" >> tmp.tmp
	echo "${SAMPLE2}_${SAMPLE3} ${n23}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE3} ${n13}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2}_${SAMPLE3} ${n123}" >> tmp.tmp
	awk '{if ($2=="") print $1"\t""0"; else print $1"\t"$2}' tmp.tmp >| ${REGION}_overlapping_peaks
	rm -f *.tmp
	Rscript ${RDIR}/make_venn_diagram_with_lines.R ${REGION}
done

### quantify number of peaks per gene
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	cut -f4  ${REGION}.peaks_present.in.at.least.2.replicates.bed | awk 'BEGIN {FS="_"} {print $1}' | sort | uniq -c  |awk 'BEGIN {print "gene.id\tno.peaks"} {FS=" "} {print $2"\t"$1}'>| ${REGION}.peaks_present.in.at.least.2.replicates_peaks.per.gene
done

### retrieve sequences at replicated peaks; make a graph of nucleotides represented at peaks; graphs produced assume that the peaks are at crosslinked sites, which may not be true
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	mkdir sequence_analysis
	cd sequence_analysis
	bedtools getfasta -s -fi /lab/solexa_page/maria/genome_files/mm10_chrMT.fa -bed ${MERGED}/${REGION}/${REGION}.peaks_present.in.at.least.2.replicates.bed |grep -v ">"| sed 's/\(.\)/\1\n/g'| sed '/^$/d' | sed 's/.*/\U&/' >| ${REGION}.peaks_present.in.at.least.2.replicates_CL.nts.only ### format is one nucleotide per line, representing a crosslinked site
	Rscript /lab/solexa_page/maria/scripts/CLIP_analysis/nucleotides_at_crosslinked_sites.R ${REGION}
done






#############
### repeat for IgG control samples; identify replicated peaks in IgG controls
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only 
	do
	PEAKS="${DIR1}/foundpeaks/demultiplexed_${IGG1}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="\t"; OFS="\t"; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12}' mm10_GRCm38_${REGION}.reg.peaks| sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12 -o distinct,max,max,max,min,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12}' | uniq >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed .

done

### for peaks that are found in multiple "regions", remove peaks based on the following hierarchy: 3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### for example, if a peak shows up in both 3' UTR and ncRNA, the hierarchical priortization removes it from ncRNA and leaves it in 3' UTR
### if there is only partial overlap in the peaks, the full peak will be removed

for REGION in intron_coding.only; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${IGG1}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in retrogenes; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done

for REGION in ncRNA.no.rRNA; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done

for REGION in cds; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;


done



for REGION in utr5_coding.only; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done

for REGION in utr3_coding.only tRNAs rRNA.only intergenic; do 
	cd ${DIR1}/foundpeaks/merged.peaks.191002/${REGION}
	cat ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done






### repeat for MM253 samples

cd ${DIR2}/foundpeaks
mkdir merged.peaks.191002


### convert peak files to bed files; merge peaks so that one peak can represent multiple transcripts; link to merged.peak directory
### filter peaks for p val <0.0001
### if no peaks meet the designated cut offs, bedtools merge will give the error message: "ERROR: Requested column 4, but database file stdin only has fields 1 - 0."
### bedtools output the strand in the wrong column so modifying bedtools output to follow bed format
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only 
	do
	mkdir "${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/"
	PEAKS="${DIR2}/foundpeaks/demultiplexed_${IGG2}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="\t"; OFS="\t"; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12}' mm10_GRCm38_${REGION}.reg.peaks| sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12 -o distinct,max,max,max,min,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12}' | uniq >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed .
	PEAKS="${DIR2}/foundpeaks/demultiplexed_${IGG3}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${PEAKS}
	awk '{FS="\t"; OFS="\t"; if ($12 < 0.001) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12}' mm10_GRCm38_${REGION}.reg.peaks| sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12 -o distinct,max,max,max,min,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12}' | uniq >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/
	ln -s ${PEAKS}/${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed .
done

### for peaks that are found in multiple "regions", remove peaks based on the following hierarchy: 3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### for example, if a peak shows up in both 3' UTR and ncRNA, the hierarchical priortization removes it from ncRNA and leaves it in 3' UTR
### if there is only partial overlap in the peaks, the full peak will be removed

for REGION in intron_coding.only; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${IGG2}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${IGG3}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in retrogenes; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${IGG3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >|  ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in ncRNA.no.rRNA; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${IGG3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in cds; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}


	bedtools subtract -header -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >|  ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${IGG3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done



for REGION in utr5_coding.only; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}

	bedtools subtract -header -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	bedtools subtract -header -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${IGG3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done

	### saving utr3_coding.only peak file with hierachry.applied in name to keep file names consistent; no utr3 peaks are removed
for REGION in utr3_coding.only tRNAs rRNA.only intergenic; do 
	cd ${DIR2}/foundpeaks/merged.peaks.191002/${REGION}
	cat ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	cat ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
done





cd ${DIR1}
cd ..
mkdir ASPeak_MM247_MM253_peaks_called_independently_hierarchy_applied
MERGED=/lab/solexa_page/maria/iCLIP_DAZL_lepts/20191001_ASPeak/ASPeak_MM247_MM253_peaks_called_independently_hierarchy_applied
cd ${MERGED}

### identify peaks that are present in at least 2 replicates using bedtools; one loop for regions where hierachy was applied; second loop for other regions
for REGION in cds intron_coding.only ncRNA.no.rRNA retrogenes utr3_coding.only utr5_coding.only; do
	cd ${MERGED}
	mkdir ${REGION}
	cd ${REGION}
	CLIP1=${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP2=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP3=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	ln -s ${CLIP1};
	ln -s ${CLIP2};
	ln -s ${CLIP3};
	bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG2}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${IGG2}_${IGG3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${IGG1}_${IGG2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG2}_${IGG3}_intersect_${REGION}.bed;
done

for REGION in intergenic rRNA.only tRNAs; do
	cd ${MERGED}
	mkdir ${REGION}
	cd ${REGION}
	CLIP1=${DIR1}/foundpeaks/merged.peaks.191002/${REGION}/${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP2=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	CLIP3=${DIR2}/foundpeaks/merged.peaks.191002/${REGION}/${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	ln -s ${CLIP1};
	ln -s ${CLIP2};
	ln -s ${CLIP3};
	bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG2}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${IGG2}_${IGG3}_intersect_${REGION}.bed;
	bedtools intersect -header -s -a ${IGG1}_${IGG2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${IGG1}_${IGG2}_${IGG3}_intersect_${REGION}.bed;
done


### create a master bed file with genomic position information only for peaks that are present in at least 2 of the 3 biological replicates
### I used the following strategy instead of bedtools because bedtools merge was doing funny things when I tried to merge intersected files
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
cd ${MERGED}/${REGION}/
cut -f1-3,6 ${IGG1}_${IGG2}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq >| tmp1.tmp
cut -f1-3,6 ${IGG1}_${IGG3}_intersect_${REGION}.bed| sort -k1,1 -k2,2n | uniq >| tmp2.tmp
cut -f1-3,6 ${IGG2}_${IGG3}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq  >| tmp3.tmp
cat tmp1.tmp tmp2.tmp >| tmp.tmp
cat tmp.tmp tmp3.tmp | sort -k1,1 -k2,2n | uniq | awk '{FS="\t"; print $1"\t"$2"\t"$3"\t""name""\t""score""\t"$4}'>|  IgG_${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed
rm -f tmp*
done


####pull out replicated peaks from each replicates peak.bed file 
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	PEAKS=${MERGED}/${REGION}/IgG_${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed;
	REP1="${DIR1}/foundpeaks/demultiplexed_${IGG1}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP1};
	bedtools intersect -header -s -a ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	REP2="${DIR2}/foundpeaks/demultiplexed_${IGG2}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP2};
	bedtools intersect -header -s -a ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	REP3="${DIR2}/foundpeaks/demultiplexed_${IGG3}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${REP3};
	bedtools intersect -header -s -a ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${PEAKS} >| ${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
done

### combine bare-bones master bed file with peak information to create a master bed file with more complete info

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}	
	FILE1="${DIR1}/foundpeaks/demultiplexed_${IGG1}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE2="${DIR2}/foundpeaks/demultiplexed_${IGG2}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE3="${DIR2}/foundpeaks/demultiplexed_${IGG3}_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	cat ${FILE1} ${FILE2} >| tmp1.tmp
	cat tmp1.tmp ${FILE3} | sort -k1,1 -k2,2n | bedtools merge -header -s -c 4,5,6,7,8,9,10,11,12 -o distinct,max,distinct,max,max,min,max,max,max -i stdin | uniq |awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' >| IgG_${REGION}.peaks_present.in.at.least.2.replicates.bed; # uniq command gets rid of repeated headers; bedtools behaving weird and output 2 strand columns
	rm -f tmp*;
done



### create a list of genes with peaks, without header

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	cut -f4 IgG_${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_"; if (NR!=1) print $1}' | sort | uniq >| IgG_${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


for REGION in retrogenes; do 
	cd ${MERGED}/${REGION}
	cut -f4 IgG_${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_exon"; if (NR!=1) print $1}' | sort | uniq >| IgG_${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


### calculate overlap in peaks among replicates; produce venn diagrams
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	area1="$(grep -v "#" 	${MERGED}/${REGION}/${IGG1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq | wc | awk '{print $1}')"
	area2="$(grep -v "#" ${MERGED}/${REGION}/${IGG2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	area3="$(grep -v "#" ${MERGED}/${REGION}/${IGG3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	n12="$(grep -v "#" ${MERGED}/${REGION}/${IGG1}_${IGG2}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n23="$(grep -v "#" ${MERGED}/${REGION}/${IGG2}_${IGG3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n13="$(grep -v "#" ${MERGED}/${REGION}/${IGG1}_${IGG3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n123="$(grep -v "#" ${MERGED}/${REGION}/${IGG1}_${IGG2}_${IGG3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	echo "Replicate No_peaks" >| tmp.tmp
	echo "${IGG1} ${area1}" >> tmp.tmp
	echo "${IGG2} ${area2}">> tmp.tmp
	echo "${IGG3} ${area3}" >> tmp.tmp
	echo "${IGG1}_${IGG2} ${n12}" >> tmp.tmp
	echo "${IGG2}_${IGG3} ${n23}" >> tmp.tmp
	echo "${IGG1}_${IGG3} ${n13}" >> tmp.tmp
	echo "${IGG1}_${IGG2}_${IGG3} ${n123}" >> tmp.tmp
	awk '{if ($2=="") print $1"\t""0"; else print $1"\t"$2}' tmp.tmp >| IgG_${REGION}_overlapping_peaks
	rm -f *.tmp
	Rscript ${RDIR}/make_venn_diagram_with_lines.R IgG_${REGION}
done

### quantify number of peaks per gene
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${MERGED}/${REGION}
	cut -f4  IgG_${REGION}.peaks_present.in.at.least.2.replicates.bed | awk 'BEGIN {FS="_"} {print $1}' | sort | uniq -c  |awk 'BEGIN {print "gene.id\tno.peaks"} {FS=" "} {print $2"\t"$1}'>| IgG_${REGION}.peaks_present.in.at.least.2.replicates_peaks.per.gene
done







#### motif analysis via MEME
### using canonical isoform for each gene
### background is the 3' UTRs of these isoforms

cd ${MERGED}/utr3_coding.only/sequence_analysis
mkdir meme_analysis
cd meme_analysis

ln -s /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10_GRCm38_canonical.transcripts.only_Ensembl_w.Refseq.id.txt
ln -s ../../utr3_coding.only.peaks_present.in.at.least.2.replicates_gene.ID.only
ln -s ../../utr3_coding.only.peaks_present.in.at.least.2.replicates.bed

grep -wFf utr3_coding.only.peaks_present.in.at.least.2.replicates_gene.ID.only mm10_GRCm38_canonical.transcripts.only_Ensembl_w.Refseq.id.txt | cut -f5 >| utr3_coding.only.peaks_canonical.isoform_Refseq.id

### extract background sequences
bsub -q 18 -K Rscript ${RDIR}/extract_3UTR_sequences.R 3UTR sleep 10 &
wait
sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' 3UTR.fa >| 3UTR.temp.fa
rm -f 3UTR.fa
mv 3UTR.temp.fa 3UTR.fa
awk 'BEGIN {FS="\t"} {print ">"$0}' utr3_coding.only.peaks_canonical.isoform_Refseq.id >| temp.id

grep -A1 -wFf temp.id 3UTR.fa |sed '/--/d' >| 3UTR_sequences_Dazl.targets_canonical.isoform.only.fa
sed 's/T/U/g' 3UTR_sequences_Dazl.targets_canonical.isoform.only.fa >| 3UTR_sequences_Dazl.targets_canonical.isoform.only_rna.fa

### expand bed file by 20 nucleotides on each side of crosslinked sites; extract as fasta sequences
	### "-q 18" send job to ubuntu cluster 18 nodes
### R script requries a list of transcripts from which I want to extract sequences
bsub -q 18 -K Rscript ${RDIR}/obtain_3UTRseq_from_bed_190220.R utr3_coding.only.peaks_present.in.at.least.2.replicates.bed utr3_coding.only.peaks_canonical.isoform_Refseq.id 20 temp sleep 10 &
wait

	### shorten header names so that meme can deal with them more easily; convert to RNA sequence
awk 'BEGIN {FS="\n"} {if (substr($1,1,1) == ">") print ">"NR; else print $0;} ' temp.fa | sed 's/T/U/g' >| utr3_coding.only.peaks_present.in.at.least.2.replicates_expanded20.rna_header.modified.fa



BKGDFASTA="3UTR_sequences_Dazl.targets_canonical.isoform.only_rna.fa"
bsub -q 18 -K /usr/local/meme/bin/fasta-get-markov -rna -m 0  ${BKGDFASTA} Dazl.bound.genes_canonical.isoform_Refseq.id_markov.model.0.order.txt sleep 10 &
wait

### motif analysis; note that consecutive crosslinked sites are merged so will only be analyzed once
### background: all expressed 3' UTRs (gene TPM>=1); if more than one isoform is expressed, most robustly expressed isoform used
DAZLBOUND="utr3_coding.only.peaks_present.in.at.least.2.replicates_expanded20.rna_header.modified.fa"
bsub -n 16 "/usr/local/meme/bin/meme ${DAZLBOUND} -rna -oc output_bkgd.canonical.isoform.markov.0.order_maxw6 -mod oops -nmotifs 6 -minw 3 -maxw 6 -maxsize 2000000 -p 16 -bfile Dazl.bound.genes_canonical.isoform_Refseq.id_markov.model.0.order.txt"

rm -f temp*




#### Identify whether GUU motif is enriched in DAZL binding sites for each region
for REGION in cds ncRNA.no.rRNA retrogenes utr5_coding.only utr3_coding.only intron_coding.only intergenic; do
	cd ${MERGED}/${REGION}/sequence_analysis
	mkdir GUU_analysis_ame
	cd GUU_analysis_ame
	ln -s ../../${REGION}.peaks_present.in.at.least.2.replicates.bed


	### expand bed file by 2 nucleotides on each side of crosslinked sites for GUU analysis; add line number to each name so each sequences has a uniqueID using NR
	awk '{FS="\t"} {if (NR==1) {print $0}  else {print $1"\t"$2-2"\t"$3+2"\t"$4"_"NR "\t"$5"\t"$6}} ' ${REGION}.peaks_present.in.at.least.2.replicates.bed >| temp_expanded_GUU.bed 


	### get fasta file
	FA="/nfs/genomes/mouse_mm10_dec_11_no_random/fasta_whole_genome/mm10.fa"
	bedtools getfasta -name -s -fi ${FA} -bed temp_expanded_GUU.bed  >| temp_expanded_GUU.fa

	## replace T with U
	sed 's/T/U/g'  temp_expanded_GUU.fa>| temp_expanded_GUU_rna.seq.fa


	### shorten header names so that meme can deal with them more easily
	awk 'BEGIN {FS="\n"} {if (substr($1,1,1) == ">") print ">"NR; else print $0;} ' temp_expanded_GUU_rna.seq.fa>| temp_expanded_GUU_rna.seq_header.modified.fa


	### using RNA library using pre-created RNA alphabet file in meme format http://meme-suite.org/doc/alphabet-format.html?man_type=web
	yes| cp  -i ~/bin/for_meme_motif_analysis/alphabet_rna.txt alphabet_rna.txt 
	### get GUU, UGUU(U/A) motif in meme format
	### using RNA library using pre-created RNA alphabet file in meme format http://meme-suite.org/doc/alphabet-format.html?man_type=web
	/usr/local/meme/bin/iupac2meme -alph alphabet_rna.txt GUU >| GUU_motif_meme_format
	
	## file for both motifs at once
	/usr/local/meme/bin/iupac2meme -alph alphabet_rna.txt GUU UGUUW >| motifs_motif_meme_format



	### create background of shuffled sequences
	/usr/local/meme/bin/fasta-shuffle-letters temp_expanded_GUU_rna.seq_header.modified.fa shuffled_control_GUU.fa

	### search for GUU, UGUU(U/A) motif using ame; set max p value to 1 so I can see non significant enrichment
	/usr/local/meme/bin/ame --verbose 1 --oc . --control shuffled_control_GUU.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 1 temp_expanded_GUU_rna.seq_header.modified.fa GUU_motif_meme_format
	mv -f ame.html ${REGION}_GUU_ame.html

	
	##make motif for GUU, UGUU(U/A)
	/usr/local/meme/bin/ceqlogo -i1 GUU_motif_meme_format -o GUU.png -f PNG

	
	#rm -f temp*

done


