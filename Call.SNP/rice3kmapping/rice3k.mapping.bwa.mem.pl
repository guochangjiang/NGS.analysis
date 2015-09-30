#!/usr/bin/perl -w
use strict;

##
## usage: $ perl rice3k.mapping.bwa.mem.pl -list listfile
##

my $splist=$ARGV[1];
open IN,"<$splist" or die "cannot open $splist\n";
while(<IN>){
        print ;
	chomp;
	my $dir=$_;
	$dir=~/^\/mnt\/hdd6\/00_Fastq\/(.+?)\/(.+?)\// or die "No match!";
	my ($type,$sample)=($1,$2);
	if (-e "/home/family/gcj/rice3k/flag/$sample.flag") {next;
	}else{
	   open OUT,">/home/family/gcj/rice3k/flag/$sample.flag";
	   close OUT;
	}
	mkdir "/home/family/gcj/rice3k/00.fastq/$sample";
	mkdir "/home/family/gcj/rice3k/mapping/$sample";
	mkdir "/home/family/gcj/rice3k/01_BAM/$type";
	mkdir "/mnt/hdd6/01_BAM/$type";
	
	print "Do $sample........";
	## trans SRA to BAM
	&SRA_mapping($type,$sample);
	##FQ to BAM
	&FQ_mapping($type,$sample);
	##merge BAM
	&Merge_BAM($sample);
	##Dedup BAM
	&Dedup($sample);	
	##Realign BAM
	&Realign($sample,$type);
	##delete unnessary files and copy result BAM to HDD6
	if (-e "/home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam"){
		system "cp -rf /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam /mnt/hdd6/01_BAM/$type/"; ## copy result bam to hdd6
        system "rm -rf /home/family/gcj/rice3k/00.fastq/$sample";  ##remove fastq
        system "rm -rf /home/family/gcj/rice3k/mapping/$sample";   ##remove mapping bam
    }
    ## calculate the DepthOfCoverage of sample
    mkdir "/home/family/gcj/rice3k/DepthOfCoverage";
    system "java -jar /home/family/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \\
			--omitDepthOutputAtEachBase \\
			-R /home/family/gcj/rice3k/ref/IRGSP-1.0_genome.fasta \\
			-T DepthOfCoverage \\
			-o /home/family/gcj/rice3k/DepthOfCoverage/$sample.IRGSP.bwa.DepthOfCoverage \\
			-I /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam
		   ";
	print "\n";
	
}
close IN;


sub SRA_mapping{
	my ($type,$sample)=@_;
	my @sra=glob "/mnt/hdd6/00_Fastq/$type/$sample/*.sra";
	foreach my $ tmp (@sra){
	     $tmp=~/(ERR\d+)/;
	     my $name=$1;
	     print "\n\t SRA mapping $sample -> $name .....";
	     next if -e "/home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sort.bam";
	     mkdir "/home/family/gcj/rice3k/00.fastq/$sample";
	     &fastq_dump($tmp,$sample);
	     my $fq1="/home/family/gcj/rice3k/00.fastq/$sample/$name"."_1.fastq";
	     my $fq2="/home/family/gcj/rice3k/00.fastq/$sample/$name"."_2.fastq";
	     &BWA_mem($fq1,$fq2,$name,$sample);
	}
}

sub FQ_mapping{
    my ($type,$sample)=@_;
    foreach my $fq1(glob "/mnt/hdd6/00_Fastq/$type/$sample/*_1.fq.gz"){
        my $fq2=$fq1;
        $fq2=~s/_1\.fq\.gz$/_2\.fq\.gz/;
        $fq1=~/00_Fastq\/$type\/$sample\/(.+)_1\.fq\.gz$/;
        my $name=$1;
        print "\n\t FQ mapping $sample -> $name .....";
        next if -e "/home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sort.bam";
        &BWA_mem($fq1,$fq2,$name,$sample);
    }
    
}
sub fastq_dump{
    my $sra=$_[0];
    my $samp=$_[1];
    system "fastq-dump --helicos --split-3 -Q 64 $sra --outdir /home/family/gcj/rice3k/00.fastq/$samp";
}

sub BWA_mem{
    my ($fq1,$fq2,$name,$sample)=@_;
    print "\n\t\t bwa mem $sample -> $name .......";
    unless (-e "/home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sam"){
        system "bwa mem -t 4 -M \\
				-R \"\@RG\tID:$sample\tLB:$sample\tPL:Illumina\tPU:$sample\tSM:$sample\"   \\
				/home/family/gcj/rice3k/ref/IRGSP-1.0_genome.fasta \\
				$fq1 $fq2   \\
				> /home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sam \\
				2> /home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.log
				";
    }
    
    system "java -Djava.io.tmpdir=/home/family/Data/tmp -jar /home/family/biosoft/picard-tools-1.114/SortSam.jar  \\
			INPUT=/home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sam   \\
			OUTPUT=/home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.sort.bam   \\
			SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT   \\
			>> /home/family/gcj/rice3k/mapping/$sample/$name.IRGSP.bwa.log 2>&1
			";
}

sub Merge_BAM{
    my $sample=$_[0];
    print "\n\t merge bam $sample ......";
    my @bam=glob "/home/family/gcj/rice3k/mapping/$sample/*.sort.bam";
    my $input;
    foreach my $bam(@bam){
        $input.="INPUT=$bam ";
    }
    unless (-e "/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.bam"){
        system "java -jar /home/family/biosoft/picard-tools-1.114/MergeSamFiles.jar   \\
				VALIDATION_STRINGENCY=LENIENT OUTPUT=/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.bam    \\
				SORT_ORDER=coordinate  \\
				$input  
				> /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.process.log 2>&1
				";
    
    }
    
}

sub Dedup{
    my $sample=$_[0];
    print "\n\t Dedup $sample ......";
    unless (-e "/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam"){
        system "java -jar /home/family/biosoft/picard-tools-1.114/MarkDuplicates.jar  \\
				MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 VALIDATION_STRINGENCY=LENIENT   \\
				INPUT=/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.bam   \\
				OUTPUT=/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam      \\
				METRICS_FILE=/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.metrics  \\
				>> /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.process.log 2>&1
				";
    }
    unless (-e "/home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam.bai"){
        system "samtools index /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam";
    }
    
}

sub Realign{
    my ($sample,$type)=@_;
    print "\n\t Realign $sample .......";
    unless (-e "/home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals"){
        system "java -jar /home/family/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \\
				-et NO_ET -K /home/family/gcj/rice3k/ref/evolution_smail.nju.edu.cn.key  \\
				-R /home/family/gcj/rice3k/ref/GATK/IRGSP-1.0_genome.fasta \\
				-T RealignerTargetCreator \\
				-nt 4  \\
				-o /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals \\
				-fixMisencodedQuals \\
				-I /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam  \\
				>> /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.process.log 2>&1
			   ";
    }
    unless (-e "/home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam"){
        system "java -jar /home/family/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar \\
				-et NO_ET -K /home/family/gcj/rice3k/ref/evolution_smail.nju.edu.cn.key \\
				-R /home/family/gcj/rice3k/ref/GATK/IRGSP-1.0_genome.fasta  \\
				-T IndelRealigner   \\
				-targetIntervals /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.intervals  \\
				-fixMisencodedQuals \\
				-o /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam    \\
				-I /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.sort.dedup.bam  
				>> /home/family/gcj/rice3k/mapping/$sample/$sample.IRGSP.bwa.process.log 2>&1
				";
    }
    
}
