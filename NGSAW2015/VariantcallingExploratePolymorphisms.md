## Variant calling and exploration of polymorphisms

### 获取数据和安装额外的模块

安装工具：

```bash

## bwa
wget -O bwa-0.7.10.tar.bz2 http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
tar xvfj bwa-0.7.10.tar.bz2
cd bwa-0.7.10
make

cp bwa /usr/local/bin

##packages
yum -y install samtools screen r-rcan-gplots python-matplotlib sysstat
```

下载数据：

```bash
git clone https://github.com/schimar/ngs2014_popGen.git

cd ngs2014_popGen/var_call2/
```

### vatiant calling

```bash

##index ref genome
bwa index ref_genome.fa

##mapping
bwa aln ref_genome.fa read_file.fq > mapped_reads.sai

##creat sam file
bwa samse ref_genome.fa mapped_reads.sai read_file.fq > mapped_reads.sam

##index ref genome
samtools faidx ref_genome.fa

##convert sam to bam
samtools view -b -S -o mapped_reads.bam mapped_reads.sam

##sort bam
samtools sort mapped_reads.bam mapped_read.sorted

##index bam
samtools index mapped_reads.sorted.bam

##查看比对结果
samtools tview mapped_reads.sorted.bam ref_genome.fa
```

### 使用Biocondctor进行variant探索

在R中输入以下代码：

```r
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("VariantAnnotation")
biocLite("SNPlocs.Hsapiens.dbSNP.20101109")
biocLite("BSgenome.Hsapiens.UCSC.hg19_1.3.1000")
```

### 质量控制

本练习的目的是比较dbSNP中call的SNP与新SNP的质量：

```r

#载入包和数据
library(VariantAnnotation)
f1 <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")

#载入样本数据，使用scanVcfHeader探索文件信息，发现info域中的VT和RSQ：
(hdr <- scanVcfHeader(fl))
info(hdr)[c("VT", "RSQ"),]

#输入数据并在其位置处峰化
(vcf <- readVcf(fl, "hg19"))
head(rowData(vcf), 3)

#SNP由MaCH/thunder call出，需将染色体编号由22变为chr22
rowData(vcf) <- renameSeqlevels(rowData(vcf), c("22"="ch22"))

#载入SNP库并查看我们的SNP是否存在其中
library(SNPlocs.Hsapiens.dbSNP.20101109)

destination <- tempfile()
pre <- FilterRules(list(isLowCoverageExomeSnp = function(x) {
grepl("LOWCOV,EXOME", x, fixed=TRUE)
}))
filt <- FilterRules(list(isSNP = function(x) info(x)$VT == "SNP"))
snpFilt <- filterVcf(fl, "hg19", destination, prefilters=pre, filters= filt)
vcf_filt <- readVcf(snpFilt, "hg19")

rowData(vcf)
rowData(vcf_filt)

#比较vcf和vcf_filt会发现10376 SNP在我们的vcf文件中，而794个在数据库中
inDbSNP <- rownames(vcf) %in% rownames(vcf_filt)
table(inDbSNP)
metrics <- data.frame(inDbSNP = inDbSNP, RSQ = info(vcf)$RSQ)

#最后可视化
library(ggplot2)
ggplot(metrics, aes(RSQ, fill=inDbSNP)) +
geom_density(alpha=0.5) +
scale_x_continuous(name="MaCH / Thunder Imputation Quality") +
scale_y_continuous(name="Density") +
theme(legend.position="top")
```


