#!/usr/bin/perl -w
use strict;

##
## usage: $ perl rice3k.mapping.DepthOfCoverage.pl -list listfile
##

my $splist=$ARGV[1];
open IN,"<$splist" or die "cannot open $splist\n";
while(<IN>){
    print ;
	chomp;
	my $dir=$_;
	$dir=~/^\/mnt\/hdd6\/00_Fastq\/(.+?)\/(.+?)\// or die "No match!";
	my ($type,$sample)=($1,$2);

	## calculate the DepthOfCoverage of sample
    mkdir "/home/family/gcj/rice3k/DepthOfCoverage";
    unless(-e "/home/family/gcj/rice3k/DepthOfCoverage/$sample.IRGSP.bwa.DepthOfCoverage*") {
    	print "Do $sample........";
    	system "java -jar /home/family/biosoft/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar --omitDepthOutputAtEachBase -R /home/family/gcj/rice3k/ref/IRGSP-1.0_genome.fasta -T DepthOfCoverage -o /home/family/gcj/rice3k/DepthOfCoverage/$sample.IRGSP.bwa.DepthOfCoverage -I /home/family/gcj/rice3k/01_BAM/$type/$sample.IRGSP.bwa.sort.dedup.realn.bam > /home/family/gcj/rice3k/DepthOfCoverage/DepthOfCoverage.log 2>&1";
		print "\n";
	}
	
}
close IN;