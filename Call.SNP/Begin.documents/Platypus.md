## Platypus
-----------------------------

参考：http://www.well.ox.ac.uk/platypus

### Platypus: A Haplotype-Based Variant Caller For Next Generation Sequence Data

#### 文献

Andy Rimmer, Hang Phan, Iain Mathieson, Zamin Iqbal, Stephen R. F. Twigg, WGS500 Consortium, Andrew O. M. Wilkie, Gil McVean, Gerton Lunter.  Integrating mapping-, assembly- and haplotype-based approaches for calling variants in clinical sequencing applications. [Nature Genetics (2014)](http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3036.html) doi:10.1038/ng.3036

#### 链接

- [最新版本下载](http://www.well.ox.ac.uk/bioinformatics/Software/Platypus-latest.tgz)
- [论坛与更新公告](https://groups.google.com/forum/#!forum/platypus-users)
- [运行示例](http://www.well.ox.ac.uk/platypus-examples)
- [常见问题FAQ](http://www.well.ox.ac.uk/platypus-faq)
- [完整文档](http://www.well.ox.ac.uk/platypus-doc)
- [GitHub](https://github.com/andyrimmer/Platypus)

#### 描述

Platypus是为高通量测序数据进行高效准确的变异检测而设计的工具。通过read的本地realignment和本地assembly，Platypus实现了高敏感性和高特异性。Platypus可以探测SNPs,MNPs,short indels, replacements和长达数kb deletions（使用assembly参数）。它广泛地应用于[whole-genome](http://www.ncbi.nlm.nih.gov/pubmed/?term=24463883)、[exon-capture](http://www.nature.com/ng/journal/v45/n1/abs/ng.2492.html)和[targeted capture](http://www.nature.com/nature/journal/v493/n7432/abs/nature11725.html)等数据。Platypus已在[Thousand Genomes](http://www.1000genomes.org/)和WGS500 projects中使用，并在[ Mainstreaming Cancer Genetics programme](http://www.mcgprogramme.com/)中用于临床测序试验。在对[Stampy](http://www.well.ox.ac.uk/project-stampy)和[BWA](http://bio-bwa.sourceforge.net/)的map数据进行分析中，Platypus表现良好。虽然尚未对其他map工具的map数据进行测试，Platypus应该也具有很好的表现。Platypus已用于人类、小鼠、老鼠和黑猩猩样品以及其他几乎所有二倍体生物体中进行变异检测。它同样用于人外显子组数据中对癌症中的体细胞突变和马赛克突变进行探测。

#### 功能

Platypus从BAM文件读取数据并输出单个VCF文件，其中包括检测出的变异列表、genotype call和所有样品的likelihood。它可以检测SNP、MNP和短的indel（小于1个read长度）和长的indel（数kb的deletion和约200bp的insertion）（使用本地assembly）。Platypus可以非常高效地处理大量的BAM文件数据和通过多BAM文件进行样品扩展。二重read标记、本地realignment和变异识别与过滤等可以由单个命令飞速进行。Platypus可以运行BAM文件形式的任何输入数据，但仅仅对Illumina数据进行过正确测试。

#### 依赖

Platypus用[Python](http://www.python.org/), [Cython](http://cython.org/)和C写出，仅需要Python2.6+和C编译器进行变异。这些工具在大多数linux和Mac OS中都已具有，Platypus的编译和运行应没有问题。

#### 建立Platypus

命令：

	$ tar -xvzf Platypus_x.x.x.tgz
	$ cd Platypus_x.x.x
	$ ./buildPlatypus.sh

> 当看到信息`Finished build Platypus`时就可以进行variant-calling了。

#### 运行

Platypus可以通过Python命令行进行运行。它需要一个或多个BAM文件和一个参考fasta文件。BAM文件需要通过[samtools](http://samtools.sourceforge.net/)进行index，fasta文件同样需要使用`samtools faidx`进行index。语法如下：
```
$ python Platypus.py callVariants --bamFiles=input.bam --refFile=ref.fasta --output=VariantCalls.vcf
```
输出一个log.txt文件和单个[VCF](http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41)文件，其中含有Platypus检测出的所有变异。log文件中的最后一行和命令结束信息都应该是`Finished variant calling`。这表示calling没有出错的完成了。另外，查看log文件中的警告和出错信息是很有用的。

#### 联系

问题反馈邮箱：platypus-users@googlegroups.com


### 运行示例

该内容包括platypus的几种不同运行情形。大多数情况下，Platypus可以在其默认设置下完美工作，而不需任何修改。要运行Platypus，首先确保所有BAM文件都已index（使用`samtools index`或其他等效兼容工具），也要使用`samtools faidx`对参考序列的fasta文件进行index。此外，还必须使用同一个参考文件进行多BAM文件的变异检测。

#### 全基因组数据的variant-calling

+ 使用Platypus分析全基因组数据是非常简单的，可通过默认设置进行分析：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fasta --output=out.vcf
```

+ 如果只在特定区域进行：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fasta --output=out.vcf --regions=chr1
```
+ 如果在多个指定区域进行：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fasta --output=out.vcf --regions=chr1,chr2,chr4
```
+ 如果在指定染色体的指定位置进行：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fasta --output=out.vcf --regions=chr1:0-1000,chr3:2200-4300
```
   此时，可以通过一个文本文件列出自定义区域信息，格式为：
```
chr1:100-20000
chr1:230000-500000
chr20:22228-99999999
...
```
   然后调用如下：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fasta --output=out.vcf --regions=FileofRegions.txt
```

#### 外显子捕获(exome capture)数据的variant-calling

使用Platypus对exome capture数据进行分析也是相当简单的，且和全基因组数据一样完美工作。默认设置已足够有效。唯一需要注意的就是Platypus会在所有数据上进行variant-calling，这意味着将进行全基因组范围的call，而不论序列数据的map位置。这是因为exome capture也能捕获到想但多的非外显子区域，并map back到基因组的不同位置，所以Platypus尝试去进行call而不管此处是否具有数据。如只要在外显子区域进行，需要指定外显子的区域信息，方法与全基因分析中的`--regions`一致。

#### 目标捕获数据(gene/exon panels etc)的variant-calling

目标基因或外显子panel或其他小区域的calling，Platypus也能很好的工作。如要关闭默认的重复read过滤（由于小区域测序导致高度重复，这与外显子组或全基因组测序不同，如过滤这些read或大大降低覆盖度），此时可以使用选项`--filterDuplicates=0`：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fa --output=out.vcf --filterDuplicates=0
```
由于目的测序具有显著偏好，所以在变异检测中谨慎对待过滤是十分必要的。一些真正的变异可能具有低于期望的等位频率（例如杂合变异仅有15%-20%的read展示），此时可能需要启用选项`--alleleBias`进行过滤。

#### Haloplex数据
如分析[Haloplex数据](http://www.genomics.agilent.com/en/Custom-NGS/HaloPlex-Custom-Kits/?cid=AG-PT-124&tabId=AG-PR-1067)，那么有一些问题需要注意。Haloplex测序数据具有大块(block)的read，其中所有的read对起始和终止位置相同，这意味着该block中的所有read都会被Platypus当作重复，所以必须Platypus的重复过滤。此外，Platypus有效去除read末端，并忽略出现在read末端的变异。修剪功能有参数`minFlank`控制。处理Haloplex数据的用法：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fa --output=out.vcf --filterDuplicates=0 --minFlank=0
```

####已知等位基因的Genotyping
如有一个变异型列表，Platypus可以进行这些等位基因的genotype：
```
python Platypus.py callVariants --bamFiles=data.bam --refFile=ref.fa --output=out.vcf --source=listOfVariants.vcf.gz --minPosterior=0 --getVariantsFromBAMs=0
```
其中，需要提供一个使用bgzip压缩并通过tabix index过的VCF文件（bgzip和tabix可从[这里](http://samtools.sourceforge.net/tabix.shtml)获取）。产生所需VCF文件的方法为：
```
bgzip listOfVariants.vcf
tabix -p vcf listOfVariants.vcf.gz
```
其中，第一行命令产生一个压缩的VCF文件(`.vcf.gz`)，第二行命令对刚刚产生的文件进行index(.vcf.gz.tbi)。这样就可以让Platypus问询VCF的特定区域了。选项`--minPosterior=0`设置所报道的变异的最小质量得分为0.这样Platypus就可以报道这些变异体的genotype了（这些数据可能在你的输出vcf中没有展示）。如不想在输出中查看reference genotype，可以移除该选项。选项`--getVariantsFromBAM=0`阻止Platypus进行正常的variant calling，而保证输入VCF中的等位基因才被报道。如向进行正常的variant calling和genotyping，可以移除该选项（此时会允许你使用外部变异体列表以扩大产生候选，对于call大的indel有用处）。

### Platypus Documentation

Platypus variant-caller的完整用户指导手册。

#### 输入选项
Platypus具有很多可用的命令行参数，可以通过`$ python Platypus.py callVariants --help`进行查看。下面的选项是一些重要的命令行参数：

+ 一般variant calling参数

选项|控制内容|默认值
----|-----|----
-h, --help|查看callVariants的帮助信息|\
--output=OUPUT,-o OUTPUT|输出VCF文件名|AllVariants.cvf
--refFile=REFFILE|index的参考fasta文件|\
--bamFiles=BAMFILES|BAM文件列表，列表以逗号分隔，或为每个BAM名单行文本文件|\
--regions=REGIONS|进行variant calling的区域，也可以逗号分隔的列表（chr:start-end,可仅给出染色体号）形式或行分隔的文本文件形式给出|所有区域（由BAM文件头获得）
--skipRegionsFile=SKIPREGIONSFILE|不进行variant calling的区域，形式与--regions已知|\
--assemble|是否以assembler来产生候选单倍型|0
--source=SOURCEFILE|进行genotyping输入的VCF文件（列表）|\
--nCPU=NCPU|运行Platypus使用的线程数|1
--logFileName=LogFile|记录信息输出文件名|log.txt
--bufferSize=BUFFERSIZE|一次读入内存的基因组区域大小，增加该值可以缩短运行时间但增加内存占用|100000
--minReads=MINREADS|最少所需reads数目|\
--maxReads=MAXREADS|窗口中最大覆盖度|\
--verbosity=VERBOSITY|记录信息的冗余度|\
--maxReadLength=RLEN|最大read长度|\
--compressReads=COMPRESSREADS|是否压缩read，设为1可减少内存占用但速度减慢|0
--maxSize=MAXSIZE|


