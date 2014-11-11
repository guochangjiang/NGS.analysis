
## samtool

### 安装
下载：http://sourceforge.net/projects/samtools/files/samtools/

安装：解压后`$ make`,后将安装目录添加到`~/.bashrc`:`export PATH=$PATH:~/biotools/samtools-1.x`

### samtools简介
samtools由一系列用于处理Sequence Alignment/Map (SAM/BAM) 格式文件的程序构成，


### 主要功能和语法
```
samtools view -bt ref_list.txt -o aln.bam aln.sam.gz

samtools sort aln.bam aln.sorted

samtools index aln.sorted.bam

samtools idxstats aln.sorted.bam

samtools view aln.sorted.bam chr2:20,100,000-20,200,000

samtools merge out.bam in1.bam in2.bam in3.bam

samtools faidx ref.fasta

samtools pileup -vcf ref.fasta aln.sorted.bam

samtools mpileup -C50 -gf ref.fasta -r chr3:1,000-2,000 in1.bam in2.bam

samtools tview aln.sorted.bam ref.fasta

bcftools index in.bcf

bcftools view in.bcf chr2:100-200 > out.vcf

bcftools view -vc in.bcf > out.vcf 2> out.afs
```

#### view
语法：
`samtools view [-bchuHS] [-t in.refList] [-o output] [-f reqFlag] [-F skipFlag] [-q minMapQ] [-l library] [-r readGroup] [-R rgFile] <in.bam>|<in.sam> [region1 [...]`



#### faidx
语法:
`samtools faidx <ref.fasta> [region1 [...]]`
为fasta格式参考序列建立索引，或者从具索引的参考序列中提取子列。若未指定区域，那么faidx将为整个文件建立索引并创建文件`<ref.fasta.fai>`.若指定
Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence. If no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk. If regions are speficified, the subsequences will be retrieved and printed to stdout in the FASTA format. The input file can be compressed in the RAZF format.


