2014.10.31


### 主要功能

```
bwa index ref.fa
bwa mem ref.fa reads.fq > aln-se.sam
bwa mem ref.fa read1.fq read2.fq > aln-pe.sam
bwa aln ref.fa short_read.fq > aln_sa.sai
bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam
bwa sampe ref.fa aln_sa1.sai aln_sa2.sai read1.fq read2.fq > aln-pe.sam
bwa bwasw ref.fa long_read.fq > aln.sam
```
### 简介
BWA用于对大参考基因组（如人类）进行小差异度序列的Mapping操作。它具有3个算法：BWA-backtrack, BWA-SW和BWA-MEM。第一种算法为100bp以下的Illumina测序设计，后两种算法是为长序列（70bp-1Mbp）设计的，且具有相似的特点，如支持长read和分拆比对，但BWA-MEM是最新的算法，一般推荐用于高质量问询，因为其更快、更精确。其实，即使处理70-100bp的Illmunina reads，MEM算法也比backtrack精确。

对于所有算法，BWA首先都需要构建参考基因组序列的FM-index，而后比对算法以不同的命令进行工作：aln/samse/sampe(backtrack)， bwasw(SW), mem(MEM).

### 安装
下载：http://sourceforge.net/projects/bio-bwa/

加压后：`$ make`，然后将安装目录添加到`~/.bashrc`:`export PATH=$PATH:~/biotools/bwa-0.x.x`,最后使用`$ source ~/.bashrc`使生效。

### 命令和参数

#### index
**语法**:`bwa index [-p prefix] [-a algoType] <in.db.fasta>`

**参数**:

-p STR	 输出数据库的名称前缀 [默认与数据库同名]

-a STR	 构建BWT index的算法[默认为is算法]：

	is	 IS linear-time algorithm.需要5.37倍于基因组数据库的内存，且不能超过2GB。IS速度中等，为默认算法，因为其简单。
	bwtsw	 该算法用于BWT-SW，如分析人类基因组时。

index可简写为`bwa index dir/ref.fasta`

#### mem
语法:`bwa mem [-aCHMpP] [-t nThreads] [-k minSeedLen] [-w bandWidth] [-d zDropoff] [-r seedSplitRatio] [-c maxOcc] [-A matchScore] [-B mmPenalty] [-O gapOpenPen] [-E gapExtPen] [-L clipPen] [-U unpairPen] [-R RGline] [-v verboseLevel] db.prefix reads.fq [mates.fq]`

用于比对70bp-1Mbp的query序列，works by seeding alignments with maximal exact matches (MEMs) and then extending seeds with the affine-gap Smith-Waterman algorithm (SW).

如缺失mates.fq文件和-p参数，那么自动视为单末端reads。如果文件mates.fq存在，那么命令将假设reads.fq中的第i个read与mates.fq中的第i个reads为一个pair。如使用-p，那么命令将假设reads.fq中的第2i个read与2i+1个read为一个pair，此时文件mates.fq被忽略。在pair-end模式下，mem命令会根据reads推断read方向和insert size分布.

mem算法进行本地比对，对于一个问询序列的不同部分可能产生多重基本匹配结果。这也是长问询序列的一个特点。而一些工具如picard markDuplicates不能进行分拆比对，此时可使用-M参数对短分拆序列的hit作为二级匹配。

参数:

-t INT 线程数[默认值1]  
-k INT 最小种子长度（短于INT的匹配被忽略.一般比对速度对此不敏感，除非INT显著偏离20）[默认值19]  
-w INT 带宽，长于INT的gap不再存在。其实，最大gap长度还受得分矩阵和hit长度的影响，而不是单单由该值决定。[默认值100]  
-d INT Off-diagonal X-dropoff (Z-dropoff).当最佳延伸得分和当前延伸得分差值大于|i-j|*A+INT时，终止延伸，其中i和j分别为query和reference的当前位置，A为匹配得分。Z-dropoff与BLAST的X-dropoff类似，只是不对gap进行罚分。Z-dropoff不仅可以避免不必要的延伸，还可以减少长且优质匹配中的劣质比对。[默认值100]  
-r FLOAT 对长于最小种子长度*FLOAT的MEM进行种子再生成，这是一个重要的绩效调谐的启发式参数。大的值产生少的种子，导致更快的比对速度和较低的准确度。[默认1.5]  
-c INT 舍弃在基因组中出现INT次以上的MEM，是一个不敏感参数。[默认10000]  
-A INT 匹配得分[默认1]  
-B INT 不匹配罚分。序列错误率近似为{0.75*exp[-log(4)*B/A]}.[默认4]  
-O INT Gap open罚分。[默认6]  
-E INT Gap extension罚分。长度为k的gap罚分为O+k*E。[1]  
-L INT 剪辑罚分。当进行SW延伸时，bwa-mem追踪达到问询末端时的最佳得分，若该得分大于最佳SW得分减去剪辑罚分，则剪辑不被采纳。SAM的AS标签报道最佳SW得分。剪辑罚分不扣除。[默认5]  
-U INT 不配对的read对的罚分。bwa-mem对未配对read对的打分为scoreRead1+scoreRead2-INT，而对配对read对的打分为scoreRead1+scoreRead2-insertPenalty。根据两个得分比较是否强制配对。[默认9]  
-p 假设第一个输入的问询文件是隔行存储的双末端fastq文件。  
-R STR 完成read group的头文件行。'\t'可以使用被转换为SAM输出中的制表符，如"@RG\tID:lane_id\tLB:sample\tPL:Illumina\tPU:sample\tSM:sample"  
-T INT 不输出得分小于INT的比对，仅影响输出。[默认30]  
-a 输出所有发现的匹配，包括单末端和未成对的双末端。这些匹配将标记为二级比对。  
-C 将fasta/q内容附加到SAM输出。该选项可用于转移read的标记信息，如barcode到SAM输出。FASTQ内容必须遵守SAM要求，否则导致错误的SAM输出。  
-H 在SAM输出中使用硬剪辑。  
-M 标记较短的分拆hit作二次匹配（为picard兼容性）。  
-v INT 控制输出的冗长水平。该选项还未能在bwa中完整支持。理论上，值0禁止一些输出到stderr，而值1则仅输出错误，值2输出警告和错误，3输出所有正常信息，4或更大值用于调试。如为4则输出不是sam格式。[3]  

#### aln
语法: `bwa aln [-n maxDiff] [-o maxGapO] [-e maxGapE] [-d nDelTail] [-i nIndelEnd] [-k maxSeedDiff] [-l seedLen] [-t nThrds] [-cRN] [-M misMsc] [-O gapOsc] [-E gapEsc] [-q trimQual] <in.db.fasta> <in.query.fq> > <out.sai>`

#### samse
语法:
`bwa samse [-n maxOcc] <in.db.fasta> <in.sai> <in.fq> > <out.sam>`

#### sampe
语法:
`bwa sampe [-a maxInsSize] [-o maxOcc] [-n maxHitPaired] [-N maxHitDis] [-P] <in.db.fasta> <in1.sai> <in2.sai> <in1.fq> <in2.fq> > <out.sam>`

#### 
语法:
`bwa bwasw [-a matchScore] [-b mmPen] [-q gapOpenPen] [-r gapExtPen] [-t nThreads] [-w bandWidth] [-T thres] [-s hspIntv] [-z zBest] [-N nHspRev] [-c thresCoef] <in.db.fasta> <in.fq> [mate.fq]`

