## VCF( Variant Call Format，变异识别格式)文件格式说明
=======================
### 参考
+ [The Variant Call Format (VCF) Version 4.2 Specfication](http://vcftools.sourceforge.net/specs.html)
+ [SAMtools GitHub](https://github.com/samtools/hts-specs)

### 1. VCF说明

VCF是一种文本文件格式，包括元信息(meta-information)行、头文件行和数据行（基因组位置相关信息）。该格式也可以包含多样本在每个位置的基因型信息。

#### 1.1 一个例子
```
##fileformat=VCFv4.2
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3
```
该示例中按次序分别显示了：1个好的SNP、1个被过滤掉的可能的SNP(质量值小于10)、具有两个等位态的位点(T为祖先位点)、1个与参考一致的单倍型位点（如无可选等位态）、1个具有两个可选等位态的微卫星（1个2bp的deletion，1个1bp的insertion）。基因型数据由三个样品给出，其中2个是phased，第3个是unphased。每个样品显示基因型质量值(GQ)、覆盖度(DP)和单倍型质量值(HQ).这3个值以逗号分隔后接数值只有phased样品才具有。其中，微卫星call位点是unphased。

#### 1.2 元信息行

文件元信息包括##字符串及其后以`key=value`对表示的内容。推荐在元信息行中使用`INFO,FILTER和FORMAT`等条目于VCF文件主体。虽然这些条目是可选的，但在元信息行中使用之可以在格式上实现完全合乎规则。

一般元信息行需要包含的内容有如下这些：

##### 1.2.1 文件格式
`fileformat`总是**必需**的，并且必须位于文件第一行，用于指明VCF格式的版本号，例如：  
`##fileformat=VCFv4.2`

##### 1.2.2 信息域格式
`INFO`域应按如下格式进行表述，其中前4个关键词是**必须**的，而后2个则推荐使用。
`##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="source",Version="version">`

信息域各个类的可能取值类型为：整型、浮点型、Flag、字符和字符串。其中，`Number`条目是一个整型数值，描述INFO域中包含的值的数目，例如：

INFO仅包括1个数值，该值为1，而若表述一对数值，那么该值为2，等等。该值还有一些特殊的字符值：

+ 如该域每个可选等位基因一个值，那么该值为A；
+ 如该域每个可能等位基因（包括参考序列）仅有一个值，那么该值为R；
+ 如该域每个可能的基因型（与FORMAT标签相关）仅有一个值，那么该值为G；
+ 如可能值的数量变化、未知或无限时，该值为"."

`Flag`类型表示INFO域不包括一个数值条目，因此`Number`应为0.`Description`值必须由双引号包括。双引号字符可以由反斜线进行转义。`Source`和`Version`值同样应该由双引号进行包括，指定注释来源（不区分大小写，如"dbsnp")和准确的版本（如"138"）。

##### 1.2.3 过滤域格式

FILTERs应该以以下格式进行描述：
```
##FILTER=<ID=ID,Description="description">
```

例如：

```
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
```

##### 1.2.4 个体格式域的格式
FORMAT域的格式为：
```
##FORMAT=<ID=ID,Number=number,Type=type,Description="description">
```

##### 1.2.5 可选等位基因域格式

描述不确定的结构变异的符号格式为：
```
##ALT=<ID=type,Description=description>
```
其中，ID值表述结构变异的类型，可以是冒号分隔的类型或亚类型列表。ID值是大小写敏感的字符串，不包含空格或尖括号。第一层次的类型必须是以下中的一个：

+ DEL：相对于参考序列是一个deletion；
+ INS：相对于参考序列是一个insertion；
+ DUP：相对于参考序列是一个拷贝数的增加；
+ INV：相对于参看序列是一个倒置；
+ CNV：拷贝数变化区（可以同时是deletion和重复duplication）。注意：当更确定的类型可以确认时，不要使用CNV。亚类型有：
  - DUP：TANDEM 串联重复；
  - DEL:ME 相对于参考序列，可移动元件(Mobile element)的deletion；
  - INS:ME 相对于参考序列，可移动元件的insertion；

另外，强烈推荐使用一些其他标签用于描述参考序列及其他数据来源信息。这些标签都是基于SAM语法的SQ域。所有便签都是可选的（参见上述示例）。

对于所有`##INFO,##Format,##FILTER和##ALT`元信息域，额外域的可以包含于其后，例如：
```
##INFO=<ID=ID,Number=number,Type=type,Description="description",Source="description",Version="128">
```
上例中，额外项`Source`和`Version`等可选项应以字符串的形式给出，即使是数值。

##### 1.2.6 比对assembly域格式

结构变异的断点比对可以以外部文件给出：`##assembly=url`.URL指明包含有锻炼比对信息的fasta文件位置。


##### 1.2.7 Contig域格式
VCF文件中涉及的contig信息的描述，格式为：
```
##contig=<ID=ctg1,URL=ftp://somewhere.org/assembly.fa,...>
```

#### 1.2.8 样品域格式
描述VCF文件中涉及的样品信息，格式为：
```
##SAMPLE=<ID=S_ID,Genomes=G1_ID;G2_ID; ...;GK_ID,Mixture=N1;N2; ...;NK,Description=S1;S2; ...;SK>
```

#### 1.2.9 系谱域pedigree格式
记录不同基因组间的家谱关系，格式为：
```
##PEDIGREE=<Name_0=G0-ID,Name_1=G1-ID,...,Name_N=GN-ID>
or a link to a database:
##pedigreeDB=<url>
```

### 1.3 header line语法
行头中必须包含8个固定的条目：
```
1. #CHROM(染色体)
2. POS(位置)
3. ID
4. REF
5. ALT
6. QUAL
7. FILTER
8. INFO
```
如果基因型(genotype)数据存在，那么在`INFO`后跟`FORMAT`条目，之后再接任意数量的样品编号。
> 注：分隔符为`Tab`.

### 1.4 数据行

#### 1.4.1 固定域
每条记录共有8个固定域，所有数据行都是以制表符为分隔符的。其中缺少的值统一以点号`.`标出。8个固定域分别为：

  + 1.CHROM - 染色体或contig编号:fasta参考序列中的序列名，其名字需连续，且无空格无冒号。（**必需**，无空格的字符串）
  + 2.POS - 位置：参考序列中的位置，起始碱基为1.同一CHROM中，位置按数字从小到大排列。同一CHROM的同一POS只能具有一条记录行。此外，端粒区域使用0或者N+1表示，其中N表示相应CHROM序列的长度。（**必需**，整型数值）
  + 3.ID - 标识号：分号分隔的唯一标识号的列表。dbSNP variant（rs数字编号），同一标识号仅能出现一次。如无可标识号，使用点号替代。（**必需**，无空格或分号分隔）
  + 4.REF - 参考序列碱基：每个碱基必需是`A,C,T,G,N`中的一个（大小写都可以），多碱基也是允许的（如微卫星位点）。其中该字符串的首个碱基对应域前面的POS显示的位置。对于简单的indels来说，REF和ALT其一应是空白的，并且二者的字符串必需包含indel事件前的碱基，除非发生在参考序列的第一个碱基，这种情况下需要包含事件后的那个碱基。在复杂的替换或其他事件中（所有等位位点字符串中具有至少一个碱基）该填补碱基是不需要的（虽然被允许）。如果ALT等位点是一个symbolic allele（尖括号ID字符串），那么填补碱基是需要的，POS信息指明多态性前一碱基的坐标。（**必需**，字符串）
  + 5.ALT - 可选碱基：从至少一个样品中call出的非参考序列等位碱基列表（以逗号分隔）。碱基字符串由ATGCN*构成（大小写不敏感）或是尖括号形式的ID字符串（<ID>）或者一个结束替代字符串。*等位位点保留以示该位点因上游的deletion而缺失。如无可选位点（均与参考位点一致），那么该位点使用缺省值点号表示。保存等位位点字符串是不需要处理VCF文件的工具的，而大小写敏感的ID字符串则需要。（字符串，无空格，逗号分隔或尖括号ID）
  + 6.QUAL - 质量值：ALT中声明位点的Phred-scaled质量值（$-10log_{10} P$）。如unknown，缺失值应指明。（数值）
  + 7.FILTER - 过滤器状态：`PASS`表示该位置通过所有过滤器（可靠）；否则该位点未通过所有过滤器，由分号分隔的代码列表指明失败的过滤器，如`q10;s50`表示该位置的质量值在10以下（错误率大于10%）、所有样品中低于50%的样品具此数据。预留的`0`不宜作为过滤器字符串。如果过滤器为被使用，该域应设为缺省值（点号）。（字符串，无空白，分号分隔）
  + 8.INFO - 附加信息：（字符串，无空白，分号分隔，允许=，值列表可用逗号分隔）。INFO域由一系列短keys及其可选值以以下格式给出：`key=data[,data]`（分号分隔），如`NS=2;DP=10;AF=0.333,0.667;AA=T;`.任意的key是允许的，但以下的这些key是预留的（可选）：
    - AA:祖先位点
    - AC:每一个ALT等位点的genotype计数（与ALT同序列出）
    - AF:每一个ALT等位点的频率（与ALT同序列出）
	- AN:在call出的genotype中所有等位点的总数
	- BQ:该位置的RMS(均方根)碱基质量
	- CIGAR:描述如何讲可选等位点比对到参考等位位点的cigar字符串
	- DB:dbSNP membership
	- DP:combined depth across samples（样品覆盖度）
	- END:该记录中描述变体的结束位置
	- H2:membership in hapmap2
	- H3:membership in hapmap3
	- MQ:RMS mapping质量
	- MQ0:覆盖到该记录的MAPQ==0的read数量
	- NS:有数据的样品数
	- SB:该位置的链偏好
	- SOMATIC:对癌基因组来说，表示该记录是一个体细胞突变
	- VALIDATED:经由后续实验验证
	- 1000G:1000Genomes中的membership
  > INFO亚域中key可以没有相应的赋值以表述group membership（组成员关系），比如`H2`表示he SNP is found in HapMap 2.不必列出一个位点不具有的属性key。

#### 1.4.2 基因型genotype域
如果基因型信息存在，那么同一类型数据必须为所有样品都呈现出来。首先，一个FORMAT域用于指出数据类型和次序（冒号分隔的字母数字字符串）。然后接每一个样品的相应数据信息（由冒号分隔）。第一个亚域必须是基因型（GT,如果存在）。此外无必需亚域。和INFO域一样，基因型域也有一些预留的常用关键词：
  + GT:基因型，有斜线/或竖线|分隔的数值（数值个数与倍体型一致）表示。如等位位点数值为0则表示为参考等位位点基因型（REF域），1是ALT域中的第一个等位位点，2是第二个，以此类推。对于二倍体call来说，可以是0/1,1|0或1/2等等。对单倍体call来说，只有一个等位位点数值出现；三倍体call则可能是0/0/1等。前后数值一致为纯合体，否则为杂合体。如果一个cal无法给出该样品再次位置的基因型信息（如为覆盖到）可用点号指明。比如./.可能表示二倍体call未覆盖到该样品。此外，斜线表示Genotype unphased，而竖线表示genotype phased。
  + DP:该样品在此位置的read深度
  + FT：样品基因型过滤器以指明该基因型是`called`。与FILTER类似，`PASS`表示通过所有过滤器，而失败的过滤器以分号列表给出，或者用点号表示未使用过滤器。这些过滤器代码值应在元信息中给出。（字符串，无空白和分号分隔）
  + 

  + GQ:conditional genotype quality, encoded as a phred quality −10log10 p(genotype call is wrong, conditioned on the site’s being variant) (Integer)
  + HQ:haplotype qualities, two comma separated phred qualities (Integers)