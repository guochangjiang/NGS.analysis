## VCF( Variant Call Format)文件格式说明
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
`fileformat`总是必需的，并且必须位于文件第一行，用于指明VCF格式的版本号，例如：  
`##fileformat=VCFv4.2`

##### 1.2.2 信息域格式
`INFO`域应按如下格式进行表述，其中前4个关键词是必须的，而后2个则推荐使用。
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

##### 1.2.6 比对域格式

结构变异的断点比对可以以外部文件给出：`##assembly=url`.URL指明包含有锻炼比对信息的fasta文件位置。


##### 1.2.7 Contig域格式

