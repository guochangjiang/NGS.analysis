#!/usr/bin/perl -w
use strict;
use File::Find::Rule;

open IN,"File_dir.hdd6.list";
while (<IN>){
    chomp;
    my ($type,$sample,$bam)=split/,/;
    if (-e "./Flag/$sample.flag"){
        next
    }else{
        open OUT,">./Flag/$sample.flag";
        close OUT;
    }
    print "Do $type $sample\n";

    system "java -jar /home/Sheep/Desktop/program/GenomeAnalysisTK-3.2-0/GenomeAnalysisTK.jar -et NO_ET -K /home/wl/Data/biosoft/evolution_smail.nju.edu.cn.key  -R /home/Sheep/Data/Mapping_of_Oryza_Sativa/Os-Nipponbare-Ref-IRGSP-1.0_GATK/IRGSP-1.0_genome.fasta -L rice.bed   -T DepthOfCoverage    -o ./New/$sample    -I $bam --omitDepthOutputAtEachBase  >> ./New/$sample.depth.log 2>&1"
 
}