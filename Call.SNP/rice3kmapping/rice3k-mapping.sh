#!/bin/sh

##
## initial mapping with bwa mem (version 0.7.9a-r786)
##
find /home/wl/Data/data/rice/00.rawdata/rice3k/01.set1 -name "*_1.fq.gz" | \
    sed 's/_1.fq.gz$//' | xargs -n 1 -P 4 -I PREFIX \
    sh -c '
        folder=`dirname PREFIX`
        sample=`basename ${folder}`
        lane_id=`basename PREFIX`
        
        echo "[`date`]: Start mapping ${sample}: ${lane_id} ... "
        
        read1=PREFIX"_1.fq.gz"
        read2=PREFIX"_2.fq.gz"
        
        ## Align reads with BWA-MEM algorithm
        bwa mem -t 4 -M -R "@RG\tID:${lane_id}\tLB:${sample}\tPL:Illumina\tPU:${sample}\tSM:${sample}" \
            /home/wl/Data/data/rice/ref/IRGSP-1.0_genome.fasta ${read1} ${read2} \
            > PREFIX.IRGSP.bwa.sam 2> PREFIX.IRGSP.bwa.log
        
        ## sort bam file
        java -Djava.io.tmpdir=/home/wl/Data/tmp -jar /home/wl/Data/biosoft/picard-tools-1.114/SortSam.jar \
            INPUT=PREFIX.IRGSP.bwa.sam \
            OUTPUT=PREFIX.IRGSP.bwa.sort.bam \
            SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
            >> PREFIX.IRGSP.bwa.log 2>&1 && \
            rm -v PREFIX.IRGSP.bwa.sam
    '

##
## process of mapping results
##
find /home/wl/Data/data/rice/00.rawdata/rice3k/01.set1 -maxdepth 1 -mindepth 1 -type d | \
    xargs -n 1 -P 4 -I FOLDER \
    sh -c '
        sample=`basename FOLDER`
        
        echo "[`date`]: Start merge bam files ${sample} ... "
        
        INPUT_BAMs=`find FOLDER -name "*.IRGSP.bwa.sort.bam" -print | \
            xargs -I BAM_FILE echo -n "INPUT=BAM_FILE "`
        
        ## merge bam files for each sample
        java -jar /home/wl/Data/biosoft/picard-tools-1.114/MergeSamFiles.jar \
            VALIDATION_STRINGENCY=LENIENT \
            OUTPUT=/home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.bam \
            SORT_ORDER=coordinate ${INPUT_BAMs} \
            > /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.process.log 2>&1
        
        echo "[`date`]: Start marking duplicates ${sample} ... "
        
        ## mark duplicates
        java -jar /home/wl/Data/biosoft/picard-tools-1.114/MarkDuplicates.jar \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT \
            INPUT=/home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.bam \
            OUTPUT=/home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam \
            METRICS_FILE=/home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.metrics \
            >> /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.process.log 2>&1 && \
            rm -v /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.bam
        
        ## index bam file
        samtools index /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam
        
        
        echo "[`date`]: Start realigning ${sample} ... "
        
        ## realignment
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
            -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
            -R /home/wl/Data/data/rice/ref/IRGSP-1.0_genome.fasta \
            -T RealignerTargetCreator -nt 4 -fixMisencodedQuals \
            -o /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.realn.intervals \
            -I /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam \
            >> /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.process.log 2>&1
        
        java -jar /home/wl/Data/biosoft/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
            -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key \
            -R /home/wl/Data/data/rice/ref/IRGSP-1.0_genome.fasta \
            -T IndelRealigner -fixMisencodedQuals \
            -targetIntervals /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.realn.intervals \
            -o /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.realn.bam \
            -I /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam \
            >> /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.process.log 2>&1 && \
            rm -v /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam \
                  /home/wl/Data/data/rice/02.assembly/00.mapped/${sample}.IRGSP.bwa.sort.dedup.bam.bai
        
        echo "[`date`]: Finished processing ${sample}"
    '
