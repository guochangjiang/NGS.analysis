#!/usr/bin/perl -w
use strict;
use File::Find::Rule;

open IN,"File_Dir.list";
mkdir "Flag";
while (<IN>){
    chomp;
    my $dir=$_;
    #print "$dir\n";
    $dir=~/^\/mnt\/hdd2\/00_Fastq\/(.+?)\/(.+?)\// or die "No match!";
    my ($type,$sample)=($1,$2);
    if (-e "./Flag/$sample.flag"){
        next
    }else{
        open OUT,">./Flag/$sample.flag";
        close OUT;
    }
    mkdir "$sample";
    print "Do $sample\n";
    chdir "$sample";
    print "Trans SRA TO BAM\n";
    &SRA_mapping($type,$sample);
    print "Trans FASTQ to BAM\n";
    &FQ_mapping($type,$sample);
    print "Merge BAM files $sample\n";
    &Merge_BAM("./",$sample);
    print "Dedup $sample\n";
    &Dedup($sample);
    print "Realign $sample\n";
    &Realign($sample,$type);
    chdir "..";
    if (-e "/mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam"){
        system "rm -rf $sample";
    }
    
}
close IN;

sub SRA_mapping{
    my ($type,$sample)=@_;
    my @sra=glob "/mnt/hdd2/00_Fastq/$type/$sample/*.sra";
    foreach my $tmp(@sra){
        #print "$tmp\n";
        $tmp=~/(ERR\d+)/;
        my $name=$1;
        next if -e "$name.IRGSP.bwa.sort.bam";
        print "        Do SRA    $name\n";
        mkdir "$name";
        chdir "$name";
        &fastq_dump($tmp);
        chdir "..";
        my $fq1="./$name/$name"."_1.fastq";
        my $fq2="./$name/$name"."_2.fastq";
        &BWA_mem($fq1,$fq2,$name,$sample);
    }
}


sub FQ_mapping{
    my ($type,$sample)=@_;
    foreach my $fq1(glob "/mnt/hdd2/00_Fastq/$type/$sample/*_1.fq.gz"){
        my $fq2=$fq1;
        $fq2=~s/_1\.fq\.gz$/_2\.fq\.gz/;
        $fq1=~/00_Fastq\/$type\/$sample\/(.+)_1\.fq\.gz$/;
        my $name=$1;
        next if -e "./$name.IRGSP.bwa.sort.bam";
        print "        Do FQ    $name\n";
        &BWA_mem($fq1,$fq2,$name,$sample);
    }
    
}

sub fastq_dump{
    my $sra=$_[0];
    system "fastq-dump --helicos --split-3 -Q 64 $sra";
}

sub BWA_mem{
    my ($fq1,$fq2,$name,$sample)=@_;
    unless (-e "./$name.IRGSP.bwa.sam"){
        system "bwa mem -t 4 -M -R \"\@RG\tID:$sample\tLB:$sample\tPL:Illumina\tPU:$sample\tSM:$sample\"   /home/Sheep/Data/Ref/IRGSP1.0/BWA/IRGSP-1.0_genome.fasta $fq1 $fq2   > ./$name.IRGSP.bwa.sam 2> ./$name.IRGSP.bwa.log";
    }
    
    system "java -Djava.io.tmpdir=/home/Sheep/Data/tmp -jar /home/Sheep/Desktop/program/picard-tools-1.114/SortSam.jar  INPUT=./$name.IRGSP.bwa.sam   OUTPUT=./$name.IRGSP.bwa.sort.bam   SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT   >> ./$name.IRGSP.bwa.log 2>&1 ";
}


sub Merge_BAM{
    my ($dir,$sample)=@_;
    my @bam = File::Find::Rule->file()
                            ->name( '*.bam')
                            ->in( "$dir" );
    my $input;
    foreach my $bam(@bam){
        $input.="INPUT=$bam ";
    }
    unless (-e "./$sample.IRGSP.bwa.sort.bam"){
        system "java -jar /home/Sheep/Desktop/program/picard-tools-1.114/MergeSamFiles.jar   VALIDATION_STRINGENCY=LENIENT OUTPUT=./$sample.IRGSP.bwa.sort.bam    SORT_ORDER=coordinate  $input  >./$sample.IRGSP.bwa.process.log 2>&1";
    
    }
    
}
sub Realign{
    my ($sample,$type)=@_;
    unless (-e "/mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals"){
        system "java -jar /home/Sheep/Desktop/program/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar -et NO_ET -K /home/Sheep/Data/Ref/IRGSP1.0/evolution_smail.nju.edu.cn.key  -R /home/Sheep/Data/Ref/IRGSP1.0/GATK/IRGSP-1.0_genome.fasta -T RealignerTargetCreator -nt 4  -o /mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals -fixMisencodedQuals -I ./$sample.IRGSP.bwa.sort.dedup.bam  >> ./$sample.IRGSP.bwa.process.log 2>&1";
    }
    unless (-e "/mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam"){
        system "java -jar /home/Sheep/Desktop/program/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar -et NO_ET -K /home/Sheep/Data/Ref/IRGSP1.0/evolution_smail.nju.edu.cn.key -R /home/Sheep/Data/Ref/IRGSP1.0/GATK/IRGSP-1.0_genome.fasta  -T IndelRealigner   -targetIntervals /mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals  -fixMisencodedQuals -o /mnt/hdd2/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam    -I ./$sample.IRGSP.bwa.sort.dedup.bam  >> ./$sample.IRGSP.bwa.process.log 2>&1"
    }
    
    
}

sub Dedup{
    my $sample=$_[0];
    unless (-e "./$sample.IRGSP.bwa.sort.dedup.bam"){
        system "java -jar /home/Sheep/Desktop/program/picard-tools-1.114/MarkDuplicates.jar  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT   INPUT=./$sample.IRGSP.bwa.sort.bam   OUTPUT=./$sample.IRGSP.bwa.sort.dedup.bam      METRICS_FILE=./$sample.IRGSP.bwa.sort.dedup.metrics  >> ./$sample.IRGSP.bwa.process.log 2>&1";
    }
    unless (-e "./$sample.IRGSP.bwa.sort.dedup.bam.bai"){
        system "samtools index ./$sample.IRGSP.bwa.sort.dedup.bam";
    }
    
}
