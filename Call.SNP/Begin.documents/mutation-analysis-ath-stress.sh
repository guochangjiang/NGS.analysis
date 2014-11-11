#!/bin/sh
#   MAGIC_mutation.sh - Mutation analysis of arabidopsis samples from MAGIC lines.
#
#   Author: Nowind
#   Created: 2014-03-27
#   Updated: 2014-06-13
#
#   References:
#   http://www.broadinstitute.org/gatk
#   http://bio-bwa.sourceforge.net/
#   http://samtools.sourceforge.net/
#   http://www.novocraft.com/main/index.php
#
#   Change logs:
#   Version 1.0.0 14/05/12: The initial version.
#   Version 1.1.0 14/05/26: Add remain samples.
#   Version 1.1.1 14/06/13: Fixed version for archive.



VERSION=1.1.1
echo ===============[$0 v$VERSION]==============
echo `date` ": start process..."
echo



######################################################################
## Mapping and Pre-processes

bwa index /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta

samtools faidx /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta

java -jar /home/wl/Data/biosoft/picard-tools-1.114/CreateSequenceDictionary.jar \
    REFERENCE=/home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
    OUTPUT=/home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.dict


find /home/wl/Data/data/arabidopsis/01.cleandata/Ara_saline/ -name "*.sra" | \
    awk '/Saline_G10-1/' | \
    xargs -n 1 -P 4 -I SRA_FILE sh -c '
        echo "dump SRA_FILE ..."
        out_path=`dirname SRA_FILE`
        out_log=`echo SRA_FILE | sed "s/.sra/.log/"`
        fastq-dump --split-files --offset 64 --helicos --defline-seq "@\$sn/\$ri" SRA_FILE \
            -O ${out_path} 2>&1 | tee ${out_log}
    '


find /home/wl/Data/data/arabidopsis/01.cleandata/Ara_saline/ -name "*_1.fastq" | \
    sed 's/_1.fastq$//' | xargs -n 1 -P 6 -I PREFIX \
    sh -c '
        sample=`basename PREFIX`
        lane_id=`basename PREFIX`
        
        echo "[`date`]: Start mapping ${sample}: ${lane_id} ... "
        
        read1=PREFIX"_1.fastq"
        read2=PREFIX"_2.fastq"
        
        ## Align reads with BWA-MEM algorithm
        bwa mem -t 3 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
            /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta ${read1} ${read2} \
            > PREFIX.TAIR10.bwa.sam 2> PREFIX.TAIR10.bwa.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/wl/Data/tmp -jar /home/wl/Data/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=PREFIX.TAIR10.bwa.sam \
            OUTPUT=PREFIX.TAIR10.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> PREFIX.TAIR10.bwa.log 2>&1 && \
            rm -v PREFIX.TAIR10.bwa.sam

        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar /home/wl/Data/biosoft/picard-tools-1.114/MarkDuplicates.jar \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=PREFIX.TAIR10.bwa.sort.bam \
            OUTPUT=PREFIX.TAIR10.bwa.sort.dedup.bam \
            METRICS_FILE=PREFIX.TAIR10.bwa.sort.dedup.metrics \
            >> PREFIX.TAIR10.bwa.log 2>&1 && \
            rm -v PREFIX.TAIR10.bwa.sort.bam
        
        ## index bam file
        samtools index PREFIX.TAIR10.bwa.sort.dedup.bam
        
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \
            -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
            -R /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T RealignerTargetCreator -nt 4 -fixMisencodedQuals \
            -o PREFIX.TAIR10.bwa.sort.dedup.realn.intervals \
            -I PREFIX.TAIR10.bwa.sort.dedup.bam \
            >> PREFIX.TAIR10.bwa.log 2>&1
        
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \
            -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
            -R /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -T IndelRealigner -fixMisencodedQuals \
            -targetIntervals PREFIX.TAIR10.bwa.sort.dedup.realn.intervals \
            -o PREFIX.TAIR10.bwa.sort.dedup.realn.bam \
            -I PREFIX.TAIR10.bwa.sort.dedup.bam \
            >> PREFIX.TAIR10.bwa.log 2>&1 && \
            rm -v PREFIX.TAIR10.bwa.sort.dedup.bam \
                  PREFIX.TAIR10.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '

######################################################################






######################################################################
## Set1: UnifiedGenotyper SNP, start from single-sample

##
## single-sample calling
##
find /home/wl/Data/data/arabidopsis/02.assembly/01.processed/02.saline_stress/ \
    -name "*.TAIR10.bwa.sort.dedup.realn.bam" -print | \
    sed 's/.bam$//' | xargs -n 1 -I PREFIX -P 3 \
    sh -c '
        echo "Start calling PREFIX.bam ..."
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \
            -R /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
            -T UnifiedGenotyper -glm SNP -nt 6 -rf BadCigar -fixMisencodedQuals -rf MappingQuality -mmq 20 \
            -o PREFIX.ug.snp.vcf -I PREFIX.bam > PREFIX.ug.snp.log 2>&1
        echo "Finished calling PREFIX.bam"
    '

bcftools merge -O v /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/samples/*.TAIR10.bwa.sort.dedup.realn.ug.snp.vcf.gz | \
    sed 's/\.\/\./0\/0/g' | bgzip -c \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.vcf.gz

vcf_process.pl --vcf /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.vcf.gz \
    --quality 30 --rare-only 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.vcf


## detect snp mutations
cat /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.vcf | \
    perl -ne 'next if(/\#/ || /(^mi)|(^ch)/); my ($chrom, $pos, $id, $ref, $alt)=(split /\s+/)[0..4];
        next unless(($ref !~ /(\w\w)|\./) && ($alt !~ /(\w\w)|\./)); print "$chrom\t",($pos-1),"\t$pos\n";' \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.bed

##
## count reads of all alleles in each strand
##
## Base Alignment Quality (BAQ) is a new concept deployed in samtools-0.1.9+.
## It aims to provide an efficient and effective way to rule out false SNPs caused by nearby INDELs.
## The default settings for mpileup is to filter reads with bitwise flag 0X704.
## So for pileup generation the following reads will not considered at all from the bam files:
##  1) 0x0400 (aka 1024 or "d") duplicate
##  2) 0x0200 (aka 512 or "f") failed QC
##  3) 0x0100 (aka 256 or "s") non primary alignment
##  4) 0x0004 (aka 4 or "u") unmapped
find /home/wl/Data/data/arabidopsis/02.assembly/01.processed/02.saline_stress/ \
    -name "*.TAIR10.bwa.sort.dedup.realn.bam" -print | sed 's/.bam$//' | \
    xargs -n 1 -P 6 -I PREFIX \
    sh -c '
        sample=`basename PREFIX | cut -d"." -f1`
        
        echo "[`date`]: Start processing ${sample} ... "
        
        samtools mpileup -d100000 -q20 -f /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
            -l /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.bed \
            PREFIX.bam | grep -vP "\t0\t" \
            > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/readcounts/${sample}.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.MQ20.mpileup
        
        java -jar /home/wl/Data/biosoft/VarScan/VarScan.v2.3.6.jar readcounts \
            /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/readcounts/${sample}.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.MQ20.mpileup \
            --min-base-qual 20 --min-coverage 1 \
            --output-file /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/readcounts/${sample}.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.readcounts
        
        echo "[`date`]: Finished processing ${sample}"
    '


for f in `find /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/readcounts/ \
    -name "*.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.readcounts"`;
do
    library=`basename $f | cut -d"." -f1`
    sample=${library}
    
    echo "${sample} ${library} ${f}"
done > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.readcounts.list


fillVcfDepth.pl --vcf /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.vcf \
    --list /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.readcounts.list \
    --type snp --minimum-vcf --update-AD \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.vcf

cat /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.vcf | \
    detect_mutations.pl -v - --max-cmp-depth 2 --min-supp-depth 4 --max-cmp-miss 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vcf


map_records.pl --query /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.csv \
    --subject /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vcf \
    -Q 0 1 -S 0 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.vs.ug.csv

echo "snapshotDirectory /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mapping_view/gr_pub" \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.igv.txt
awk '!/#/ {print "goto "$1":"($2-50)"-"($2+50); print "snapshot "$3"_"$1"_"$2".png";}' \
    /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.csv \
    >> /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.igv.txt


map_records.pl --query /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vcf \
    --subject /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mutation.ath.gr.csv \
    -Q 0 1 -S 0 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vs.gr.csv

echo "snapshotDirectory /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/mapping_view/ug_fp" \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vs.gr.igv.txt
awk '$8 ~ /NMISS=0;SMISS=NA;FPD=0;FPFQ=0/ && $11 == "N/A" {print "goto "$1":"($2-50)"-"($2+50); print "snapshot "$3"_"$1"_"$2".png";}' \
    /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vs.gr.csv \
    >> /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.snp.rare.Q20.mut.vs.gr.igv.txt




##
## generate alignments for manually inspections
##
for record in `perl -ne 'next if (/#/); my ($chrom, $pos, $sample, $info) = (split /\t/)[0,1,2,7];
    $info =~ /MA=(\w+)/; print "$sample;$chrom:$pos#$1\n";' \
    /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/Samples_MAGIC.TAIR10.bwa.stampy.realn.ug.mut.snp.vcf`;
do
    spec_sample=${record/;*}
    mutation=${record/*;}
    mut_base=${mutation/*\#}
    mutation=${mutation/\#*}
    chrom=${mutation/:*}
    mut_pos=${mutation/*:}

    start_pos=`echo "${mut_pos}-100" | bc`
    end_pos=`echo "${mut_pos}+100" | bc`

    echo "${spec_sample} ${chrom} ${mut_pos} ${mut_base}"
    
    if [[ -n /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample} ]]; then
        mkdir -pv /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}
    fi
    
    echo -e "${chrom} ${start_pos} ${end_pos}\n${chrom} ${start_pos} ${mut_pos}\n${chrom} ${mut_pos} ${end_pos}\n" | \
        query_fasta_seqs.pl -i - \
        -f /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
        -L 0 1 2 -s 1 2 \
        > /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}/${chrom}_${mut_pos}.fa
    
    for bam_file in `ls /home/wl/Data/data/arabidopsis/02.assembly/01.processed/00.by_lib/00.normal_samples/*.TAIR10.bwa.stampy.realn.bam`;
    do
        sample=`basename ${bam_file} | cut -d"." -f1`
        samtools view ${bam_file} ${chrom}:${mut_pos}-${mut_pos} | \
            awk -v pos=${mut_pos} -v base=${mut_base} -v name=${sample} 'BEGIN {OFS = FS = "\t"}; {n=split($10,a,"");
            if(a[(pos-$4)+1] == base){ print ">"name"|"$1"\n"$10}}' \
            >> /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}/${chrom}_${mut_pos}.fa
    done

    reference_align.pl -i /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}/${chrom}_${mut_pos}.fa \
        > /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}/${chrom}_${mut_pos}.aln.fas && \
        rm -v /home/wl/Data/data/arabidopsis/03.analysis/06.mutations/distribution/snp/align/${spec_sample}/${chrom}_${mut_pos}.fa
done
######################################################################




######################################################################
## Set2: UnifiedGenotyper multi-sample

##
## multi-sample calling
##
BAM_FILES=`find /home/wl/Data/data/arabidopsis/02.assembly/01.processed/02.saline_stress/ \
    -name "*.TAIR10.bwa.sort.dedup.realn.bam" -print | xargs -I BAM echo -n "-I BAM "`
java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \
    -R /home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
    -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
    -T UnifiedGenotyper -glm BOTH -nt 12 -rf BadCigar -fixMisencodedQuals -rf MappingQuality -mmq 20 \
    -o /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.vcf \
    ${BAM_FILES} 2>&1 | \
    tee /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.ug.log


## screen out rare variants which have a high possibility of been a candidate mutation
vcf_process.pl --vcf /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.vcf \
    --quality 30 --rare-only 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.rare.vcf
######################################################################






######################################################################
## Set3: Platypus multi-sample

##
## multi-sample calling
##
BAM_FILES=`find /home/wl/Data/data/arabidopsis/02.assembly/01.processed/02.saline_stress/ \
    -name "*.TAIR10.bwa.sort.dedup.realn.bam" -print | xargs -I BAM echo -n "BAM," | sed 's/,$//'`
Platypus.py callVariants --refFile=/home/wl/Data/data/arabidopsis/ref/TAIR10_genome_release/TAIR10_chr_all.fasta \
    --bamFiles=${BAM_FILES} \
    --nCPU=12 --assemble=1 --verbosity=1 --mergeClusteredVariants=1 \
    --output=/home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.vcf \
    --logFileName=/home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.log

## screen out rare variants which have a high possibility of been a candidate mutation
vcf_process.pl --vcf /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/variants/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.vcf \
    --quality 30 --rare-only 1 \
    > /home/wl/Data/data/arabidopsis/03.analysis/09.saline_stress/mutation/Samples_saline.TAIR10.bwa.sort.dedup.realn.Platypus.rare.vcf
######################################################################






echo `date` ':all processes completed!'
echo

echo ---------------------------------------

exit 0




