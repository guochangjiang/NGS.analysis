## A complete de novo assembly and annotation protocol for mRNASeq

###下载样品数据

样品数据来自于研究[A quantitative reference transcriptome for Nematostella vectensis early embryonic development: a pipeline for de novo assembly in emerging model systems](http://www.evodevojournal.com/content/4/1/16)：

```bash
curl -O http://athyra.idyll.org/~t/mrnaseq-subset.tar
tar xvf mrnaseq-subset.tar
```

### 方法

所需工具：

https://khmer-protocols.readthedocs.org/，可用http://khmer-protocols.readthedocs.org/en/ngs2014/替代。

**流程**：

1. Quality trimming
2. Applying digital normalization
3. Running the actual assembly
4. BLASTing your assembly

详细内容（http://khmer-protocols.readthedocs.org/en/ngs2014/index.html#mrnaseq-assembly-the-eel-pond-protocol）

#### 1. 质量修剪和序列过滤

+ 安装kmer

```bash
cd ${HOME}
mkdir -p projects/eelpond
python2.7 -m virtualenv projects/eelpond/env
source ${HOME}/projects/eelpond/env/bin/activate
mkdir -p src
cd src
git clone --branch v1.3 https://github.com/ged-lab/khmer.git
cd khmer
make install
```

其中，选项`virtualenv`使我们可以在没有root权限条件下安装python软件。

+ 下载数据

[0Hour_ATCACG_L002_R1_001.fastq.gz](https://darchive.mblwhoilibrary.org/bitstream/handle/1912/5613/0Hour_ATCACG_L002_R1_001.fastq.gz?sequence=2&isAllowed=y)

为了节省时间，先不要使用整个数据进行运行，可以先提取部分数据：

```bash
cd raw
ln -fs ${HOME}/data/nemo/*.fastq.gz .
mkdir -p extract
for file in raw/*.fastq.gz
do
    gunzip -c ${file} | head -400000 | gzip \
        > extract/${file%%.fastq.gz}.extract.fastq.gz
done
```

+ 使用FastQC评估质量

```bash
fastqc --threads 1 *.fastq.gz --outdir=${HOME}/fastqc
```

+ 找到正确的Illumina adapters

需要确认是否采用了正确的adaptor，如果不是正确的adaptor，那么read将不会被修剪：

```bash
wget https://sources.debian.net/data/main/t/trimmomatic/0.32+dfsg-2/adapters/TruSeq3-PE.fa
```

+ 对每对文件进行adapter修剪：

对以下两个paired-end文件进行处理：

```
24HourB_GCCAAT_L002_R1_001.fastq.gz
24HourB_GCCAAT_L002_R2_001.fastq.gz
```

对二者分别简写为<R1 file>和<R2 file>, 且需要指定一个特定的<sample name>。

```bash
# run trimmomatic
trimmomatic PE <R1 FILE> <R2 FILE> s1_pe s1_se s2_pe s2_se \
    ILLUMINACLIP:${HOME}/projects/eelpond/TruSeq3-PE.fa:2:30:10

# interleave the remaining paired-end files
interleave-reads.py s1_pe s2_pe | gzip -9c \
    > ../trimmed/<SAMPLE NAME>.pe.fq.gz

# combine the single-ended files
cat s1_se s2_se | gzip -9c > ../trimmed/<SAMPLE NAME>.se.fq.gz

# clear the temporary files
rm *

# make it hard to delete the files you just created
cd ../trimmed
chmod u-w *
```

+ 通过脚本实现自动化

```bash
python ${HOME}/src/khmer/sandbox/write-trimmomatic.py > trim.sh
more trim.sh
bash trim.sh
```

+ 质量修剪

对于所有的`.pe.fq.gz`和`.se.fq.gz`文件，运行：

```bash
gunzip -c <filename> | fastq_quality_filter -Q33 -q 30 -p 50 | gzip -9c \
> <filename>.qc.fq.gz
````

此外，可以将命令`fastx_trimmer -Q33 -l 70 |`放入其中实现固定长度read的输出。另外，如果移除`-Q33`选项，那么可能导致质量值问题。

+ 自动化

```bash
for file in *
do
     echo working with ${file}
     newfile=${file%%.fq.gz}.qc.fq.gz
     gunzip -c ${file} | fastq_quality_filter -Q33 -q 30 -p 50 | gzip -9c \
         > ../filtered/${newfile}
done
```

+ 从间隔文件中提取paired-ends

```bash
mkdir filtered-pairs
cd filtered-pairs
for file in ../filtered/*.pe.qc.fq.gz
do
   extract-paired-reads.py ${file}
done
```

+ 结束

最后得到的文件列表：

```
raw/24HourB_GCCAAT_L002_R1_001.fastq.gz                   - the original data
raw/24HourB_GCCAAT_L002_R2_001.fastq.gz
trimmed/24HourB_GCCAAT_L002_R1_001.pe.fq.gz               - adapter trimmed pe
trimmed/24HourB_GCCAAT_L002_R1_001.se.fq.gz               - adapter trimmed orphans
filtered/24HourB_GCCAAT_L002_R1_001.pe.qc.fq.gz           - FASTX filtered
filtered/24HourB_GCCAAT_L002_R1_001.se.qc.fq.gz           - FASTX filtered orphans
filtered-pairs/24HourB_GCCAAT_L002_R1_001.pe.qc.fq.gz.pe  - FASTX filtered PE
filtered-pairs/24HourB_GCCAAT_L002_R1_001.pe.qc.fq.gz.se  - FASTX filtered SE
```

删除无用文件：

```bash
#Well, first, you can get rid of the original data. You already have it on a disk somewhere, right?

rm raw/*
rmdir raw

#Next, you can get rid of the trimmed files, since you only want the QC files. So

rm -f trimmed/*
rmdir trimmed

#And, finally, you can toss the filtered files, because you’ve turned those into *.pe and *.se files:

rm filtered/*
rmdir filtered
```

最后剩余文件：

```
filtered-pairs/24HourB_GCCAAT_L002_R1_001.pe.qc.fq.gz.pe   - FASTX filtered PE
filtered-pairs/24HourB_GCCAAT_L002_R1_001.pe.qc.fq.gz.se   - FASTX filtered SE
filtered-pairs/24HourB_GCCAAT_L002_R1_001.se.qc.fq.gz      - FASTX filtered orphans
```

注意：

1. 文件名要方便显示其信息与历史，即使丑陋也没关系
2. 每一步都可以建立一个新的目录，方便在后续步骤中删除无用文件

重命名文件：

```bash
## PE
for file in filtered-pairs/*.pe
do
   newfile=${file%%.pe.qc.fq.gz.pe}.pe.qc.fq
   mv $file $newfile
   gzip $newfile
done

##SE
for file in filtered-pairs/*.se
do
  otherfile=${file%%.pe.qc.fq.gz.se}.se.qc.fq.gz # the orphans
  gunzip -c ${otherfile} > combine
  cat ${file} >> combine
  gzip -c combine > ${otherfile} # now all the single reads together
  rm ${file} combine
done

##make the end product files read-only
chmod u-w filtered-pairs/*
```

#### 数字标准化应用

+ 连接数据

```bash
mkdir -p ${HOME}/filtered-pairs
ln -fs ${HOME}/data/*.qc.fq.gz ${HOME}/filtered-pairs/
```

+ 运行数字标准化

```bash

##paired-end reads
cd ${HOME}
mkdir diginorm
cd diginorm
normalize-by-median.py --paired -ksize 20 --cutoff 20 -n_tables 4 \
  --min-tablesize 3e8 --savetable normC20k20.ct \
  ../filtered-pairs/*.pe.qc.fq.gz

##single-end reads
normalize-by-median.py --cutoff 20 --loadtable normC20k20.ct \
  --savetable normC20k20.ct ../filtered-pairs/*.se.qc.fq.gz
```

其中，`--paired`表示处理PE数据，保证所有reads不会落单；而`--n_tables`和`--min_tablesize`参数则指定多少内存可使用，应比计算机内存小。例如，对于转录组来说，最好为~60GB左右，所以可以采用`--n_tables 4 --min_tablesize 15e9`; 而对于数百百万的reads，16GB是足够的。（详见http://khmer.readthedocs.org/en/latest/choosing-hash-sizes.html）

+ 剪掉可能错误的k-mers

遍历所有reads并剪掉高覆盖reads中的低丰度部分：

```bash
mkdir abundfilt
cd abundfilt
filter-abund.py --variable-coverage ../diginorm/normC20k20.ct \
  --threads ${THREADS:-1} ../diginorm/*.keep
```

这一步会使一部分reads落单，但是没关系，因为它们的partner read不好。

+ 重命名文件

```bash
# break out the orphaned and still-paired reads:
mkdir digiresult
cd digiresult
for file in ../abundfilt/*.pe.*.abundfilt
do
   extract-paired-reads.py ${file}
done

# combine the orphaned reads into a single file:
for file in *.se.qc.fq.gz.keep.abundfilt
do
   pe_orphans=${file%%.se.qc.fq.gz.keep.abundfilt}.pe.qc.fq.gz.keep.abundfilt.se
   newfile=${file%%.se.qc.fq.gz.keep.abundfilt}.se.qc.keep.abundfilt.fq.gz
   cat ${file} ../digiresult/${pe_orphans} | gzip -c > ../digiresult/${newfile}
   rm ${pe_orphans}
done

# rename the remaining PE reads & compress those files:
for file in *.abundfilt.pe
do
   newfile=${file%%.fq.gz.keep.abundfilt.pe}.keep.abundfilt.fq
   mv ${file} ${newfile}
   gzip ${newfile}
done
```

剩余文件：

```
filtered-reads/6Hour_CGATGT_L002_R1_005.pe.qc.fq.gz
filtered-reads/6Hour_CGATGT_L002_R1_005.se.qc.fq.gz
diginorm/6Hour_CGATGT_L002_R1_005.pe.qc.fq.gz.keep
diginorm/6Hour_CGATGT_L002_R1_005.se.qc.fq.gz.keep
diginorm/normC20k20.ct
abundfilt/6Hour_CGATGT_L002_R1_005.pe.qc.fq.gz.keep.abundfilt
abundfilt/6Hour_CGATGT_L002_R1_005.se.qc.fq.gz.keep.abundfilt
digiresult/6Hour_CGATGT_L002_R1_005.pe.qc.keep.abundfilt.fq.gz
digiresult/6Hour_CGATGT_L002_R1_005.se.qc.keep.abundfilt.fq.gz
```

#### 3. 进行正式的assembly

下面这些都会在screen中运行，至少有15GB内存

+ 安装Trinity

```bash
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/v2.0.4.tar.gz \
  -O trinity.tar.gz
tar xzf trinity.tar.gz
cd trinityrnaseq*/
make |& tee trinity-build.log
```

+ 创建进行assemble的文件

对于PE数据，Trinity需要两个文件: left & right, 它们可以含有落单reads。所以可以将间隔的pe文件拆分为两个，并将se序列添加进去：

```bash
for file in *.pe.qc.keep.abundfilt.fq.gz
do
   split-paired-reads.py ${file}
done

cat *.1 > left.fq
cat *.2 > right.fq

gunzip -c *.se.qc.keep.abundfilt.fq.gz >> left.fq
```

+ 使用Trinity进行assembly：

```bash
Trinity --left left.fq --right right.fq --seqTypr fq --max_memory 14G --CPU 1
```

最终经过长时间的运行结束后，一个组装好的转录组就会产生`Trinity.fasta`。

#### 4. 组装数据的blast

安装blastkit：

```bash
git clone https://github.com/ctb/blastkit.git -b ec2
cd blastkit/www
ln -fs $PWD /var/www/blastkit

mkdir files
chmod a+rxwt files

python check.py
```

添加数据：

```bash
cp Trinity.fasta /blastkie/db/db.fa

#或者使用已有数据
curl -O https://s3.amazonaws.com/public.ged.msu.edu/trinity-nematostella-raw.fa.gz
gunzip trinity-nematostella-raw.fa.gz
mv trinity-nematostella-raw.fa db/db.fa
```

格式化数据：

```bash
formatdb -i db/db.fa -o T -p F
python index-db.py db/db.fa
```

运行blastkit：
使用以下序列
```
CAGCCTTTAGAAGGAAACAGTGGCAATATATAATTCTAGATGAAGCTCAGAATATCAAAA
ATTTTAAAAGTCAAAGGTGGCAGTTGCTGTTGAATTTTTCAAGTCAGAGGAGACTTTTGT
TGACTGGAACACCTTTGCAGAACAATTTGATGGAGCTGTGGTCGCTTATGCATTTCCTCA
TGCCATCAATGTTTGCTTCTCATAAAGATTTTAGGGAGTGGTTTTCTAACCCTGTTACAG
GGATGATTGAAGGGAATTCAG
```

它应该会在你的组装中匹配到一些数据。

#### 5. 构建转录家族

建立个人的二进制程序文件夹：

```bash
mkdir -p ${HOME}/bin
export PATH=${PATH}:${HOME}/bin
echo 'export PATH=${PATH}:${HOME}/bin' >> ${HOME}/.bashrc
```

拷贝数据：

```bash
gzip -c trinity_out_dir/Trinity.fasta > trinity-nematostella-raw.fa.gz
```

运行khmer分隔：

Partioning运行基于de Bruijin graph的成簇算法可以通过过渡序列重叠将你的转录本簇化。也就是，将转录本分组为转录本家族（根据shared序列）：

```bash
mkdir partitions
cd partitions
do-partition.py -x 1e9 -N 4 --threads 1 nema \
  ../trinity-Nematostella-raw.fa.gz
```

过10多分钟后，输出的文件（以.part结尾）包含有分区分配信息。可以分组和重命名这些序列：

```bash
rename-with-partitions.py nema trinity-nematostella-raw.fa.gz.part
mv trinit-nematostella-raw.fa.gz.part.renamed.fasta.gz \
  trinity-nematostella.renamed.fa.gz
```

查看重命名序列

```bash
gunzip -c trinity-nematostella.renamed.fa.gz | head
```

每条序列大概具有这样的名字：

```
>nema.id1.tr16001 1_of_1_in_tr16001 len=261 id=1 tr=16001
```

解释：

+ nema	- 重命名时设定的前缀
+ idN	- 序列的独有id
+ trN	- 是转录本家族
+ 1_0f_1_in_tr16001	- 表示该转录本家族仅有一个成员
+ len	- 序列长度

#### 6. 注释转录本家族

生成所需文件（上节）或下载之：

```bash
curl -O http://public.ged.msu.edu.s3.amazonaws.com/trinity-nematostella.renamed.fa.gz
```

由于BLASTs将花费较长时间，所以可以使用canned BLASTs：

```bash
cd ${HOME}/annotation
curl -O http://public.ged.msu.edu.s3.amazonaws.com/nema.x.mouse.gz
curl -O http://public.ged.msu.edu.s3.amazonaws.com/mouse.x.nema.gz
gunzip nema.x.mouse.gz
gunzip mouse.x.nema.gz
```

针对小鼠的初步注释：

```bash
#假定同源性，做BLASTs和相互最佳匹配分析，首先解压文件
gunzip -c ../partitions/trinity-nematostella.renamed.fa.gz \
  trinity-nematostella.renamed.fa

#获取最新mouse RefSeq
for file in mouse.1.protein.faa.gz mouse.2.protein.faa.gz
do
     curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/${file}
done
gunzip mouse.[123].protein.faa.gz
cat mouse.[123].protein.faa > mouse.protein.faa

#格式化blast库
formatdb -i mouse.protein.faa -o T -p T
formatdb -i trinity-nematostella.renamed.fa -o T -p F

#在screen中进行双向blast
blastx -db mouse.protein.faa -query trinity-nematostella.renamed.fa \
   -evalue 1e-3 -num_threads 8 -num_descriptions 4 -num_alignments 4 \
   -out nema.x.mouse
   tblastn -db trinity-nematostella.renamed.fa -query mouse.protein.faa \
   -evalue 1e-3 -num_threads 8 -num_descriptions 4 -num_alignments 4 \
   -out mouse.x.nema
```

指定序列名称：根据相互最佳匹配计算假定的homology (best BLAST hit) and orthology (reciprocal best hits):

```bash
make-uni-best-hits.py nema.x.moutse nema.x.mouse.homol
make-reciprocal-best-hits.py nema.x.mouse mouse.x.nema.mouse.ortho

#准备mouse信息中的一些
make-namedb.py mouse.protein.faa mouse.namedb
python -m screed.fadbm mouse.protein.faa

#最后注释序列
annotate-seqs.py trinity-nematostella.renamed.fa nema.x.mouse.ortho \
  nema.x.mouse.homol
```

最后回看到：

```
207533 sequences total
10471 annotated / ortho
95726 annotated / homol
17215 annotated / tr
123412 total annotated
```

> 如果其中任何一个数为0，那么需要重新做BLAST。

最后生成文件`trinity-nematostella.renamed.fa.annot`，其中的序列似下：

```
>nematostella.id1.tr115222 h=43% => suppressor of tumorigenicity 7 protein isoform 2 [Mus musculus] 1_of_7_in_tr115222 len=1635 id=1 tr=115222 1_of_7_in_tr115222 len=1635 id=1 tr=115222
```

建议重命名改为该文件为`nematostella.fa`,用其进行4中的blast：

```bash
cp trinity-nematostella.renamed.fa.anoot Nematostella.fa
```

注释过程还会产生2个CSV文件，小的为`trinity-nematostella.renamed.fa.annot.csv`，包含连接到同源信息的序列名字，而大的那个`trinity-nematostella.renamed.fa.annot.large.csv`则包含所有信息。

#### 7. 表达分析

首先，获取分割和重命名数据

```bash
curl -O http://athyra.idyll.org/~t/trinity-nematostella.renamed.fa.gz
gunzip -c trinity-nematostella.renamed.fa.gz > nematostella.fa
```

+ 安装rsem

使用[RSEM](http://deweylab.biostat.wisc.edu/rsem/)进行表达分析，使用[EBSeq](http://www.biostat.wisc.edu/~kendzior/EBSEQ/)进行差异表达分析，安装方法如下：

```bash
curl -O http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.8.tar.gz
tar xzf rsem-1.2.8.tar.gz
cd rsem-1.2.8
make
cd EBSeq
make

#RSEM需要bowtie
curl -O -L http://sourceforge.net/projects/bowtie-bio/files/bowtie/0.12.7/bowtie-0.12.7-linux-x86_64.zip
unzip bowtie-0.12.7-linux-x86_64.zip
cd bowtie-0.12.7
cp bowtie bowtie-build bowtie-inspect /usr/local/bin
```

准备reference：

```bash
mkdir rsem
cd rsem

ln -fs ../nematostella.fa .

#make a transcript-to-gene-map file
python make-transcript-to-gene-map-file.py nematostella.fa nematostella.fa.tr_to_genes

#ask rsen准备reference
rsem-prepare-reference --transcript-to-gene-map nematostella.fa.tr_to_genes nematostella.fa nema
```

找到并列出reads

```bash
#找到QC reads，链接过来
ln -fs /data/*.pe.qc.fq.gz
ls -l *.pe.qc.fq.gz > list.txt
```

运行RSEM

对于`list.txt`中的每个文件运行RSEM(比较耗时间，在screen中进行)：

```bash
n=1
for filename in $(cat list.txt)
do
    echo mapping $filename
    gunzip -c $filename > ${n}.fq
    /usr/local/share/khmer/scripts/split-paired-reads.py ${n}.fq
    rsem-calculate-expression --paired-end ${n}.fq.1 ${n}.fq.2 nema -p 4 ${n}.fq
    rm ${n}.fq ${n}.fq.[12] ${n}.fq.transcript.bam ${n}.fq.transcript.sorted.bam
    n=$(($n + 1))
done
```

聚集结果

```bash
rsem-generate-data-matrix [0-9].fq.genes.results 10.fq.genes.results > 0-vs-6-hour.matrix

#也可指定顺序
rsem-generate-data-matrix 1.fq.genes.results 3.fq.genes.results > results.matrix
```

#### 8. 差异表达分析（EMSeq）


```bash
gunzip -c /data/nematostella.fa.gz > ./nematostella.fa
mkdir ebseq
cd ebseq

rsem-generate-data-matrix /data/[0-9].fq.genes.results /data/10.fq.genes.results > 0-vs-6-hour.matrix

mkdir ebseq
cd ebseq
cp ../rsem/0-vs-6-hour.matrix .

rsem-run-ebseq 0-vs-6-hour.matrix 5,5 0-vs-6-hour.changed
#.matrix file contains 2 conditions, each with 5 replicates; if you had two replicates, you would call rsem-run-ebseq with 2,2.

#提取差异表达基因并结合于注释信息
python /usr/local/share/eel-pond/extract-and-annotate-changed.py 0-vs-6-hour.changed /mnt/nematostella.fa 0-vs-6-hour.changed.csv

#可视化
python /usr/local/share/eel-pond/plot-expression.py 0-vs-6-hour.matrix 5,5 0-vs-6-hour.changed.csv
```

![](http://khmer-protocols.readthedocs.org/en/ngs2014/_images/0-vs-6-hour.matrix.png)

### The Kalamazoo Metagenome Assembly protocol

http://khmer-protocols.readthedocs.org/en/ngs2014/metagenomics/index.html#the-kalamazoo-metagenome-assembly-protocol


