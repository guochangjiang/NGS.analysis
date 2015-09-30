## 


### 亚马逊EC2入门（Ubuntu)

+ 切换管理员权限
``` 
sudo bash
```
+ 从网络获取文件
	+ wget
	+ curl -O xxx.com/file
+ 解压缩文件
	- tar -xvzf dropbox.tar.gz
	- gunzip

### BLAST使用
#### 安装前的准备工作
依赖软件包的安装：screen, git, curl, gcc, make, g++, python-dev, unzip, default-jre, pkg-config, libncurses5-dev, r-bash-core, r-cran-gplots, python-matplotlib, sysstat
```
apt-get update
apt-get -y install screen git curl gcc make g++ python-dev unzip \
        default-jre pkg-config libncurses5-dev r-base-core \
        r-cran-gplots python-matplotlib sysstat
```

#### 安装blast
```
curl -O ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.26/blast-2.2.26-x64-linux.tar.gz
tar xzf blast-2.2.26-x64-linux.tar.gz
cp blast-2.2.26/bin/* /usr/local/bin
cp -r blast-2.2.26/data /usr/local/blast-data
```

#### 使用blast
```
curl -O ftp://ftp.ncbi.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.protein.faa.gz
curl -O ftp://ftp.ncbi.nih.gov/refseq/D_rerio/mRNA_Prot/zebrafish.protein.faa.gz
gunzip *.faa.gz   #解压缩的同时删除压缩文件
head -11 mouse.protein.faa > mm-first.fa  #生成query序列文件
formatdb -i zebrafish.protein.faa -o T -p T  #对目标库文件进行建库
blastall -i mm-first.fa -d zebrafish.protein.faa -p blastp -b 2 -v 2 -e 1e-6 -o out.txt  #其中，-b 2 -v 2表示仅显示两条匹配结果，而-e 1e-6表示仅显示evalue值小于1e-6的匹配结果。
```

#### 将blast结果转换成CSV格式
+ 工具：ngs-scripts/blast/blast-to-csv.py
	使用方法：
```
git clone https://github.com/ngs-docs/ngs-scripts.git /ngs-scripts
python /ngs-scripts/blast/blast-to-csv.py out.txt > out.csv
```

### 短测序reads的质量控制

#### 软件安装

+ Trimmomatic
```
curl http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.32.zip -o /Trimmomatic-0.32.zip
        unzip -o /Trimmomatic-0.32.zip -d /Trimmomatic
        ln -fs /Trimmomatic/trimmomatic-0.32.jar ~/bin
```

+ FastQC
```
curl http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.2.zip -o /fastqc_v0.11.2.zip
        unzip -o /fastqc_v0.11.2.zip -d /FastQC
        chmod +x /FastQC/fastqc
        ln -fs /FastQC/fastqc ~/bin/fastqc
```

#### 文件准备
```
curl apollo.huck.psu.edu/data/SRR.tar.gz -o SRR.tar.gz
tar xzvf SRR.tar.gz
```

[Fastq文件格式说明](http://en.wikipedia.org/wiki/FASTQ_format)

#### 运行fastqc
```
$ fastqc *.fastqc -o /outdir
```

####  模式匹配
在文件中查找目标字符串并高亮显示
```
grep STRingxxx target.file --color=always
```

#### 修剪Trimming
使用软件Trimmomatic对read进行修剪，比如去除Illumina的adapter序列：
```
ln -s ~/Trimmomatic-0.32/adapters/TruSeq3-SE.fa
alias trim='java -jar ~/Trimmomatic-0.32/trimmomatic-0.32.jar' #赋予一个长命令以一个简写
trim SE SRR447649_1.fastq good.fq SLIDINGWINDOW:4:25 MINLEN:36 #只通过质量值进行修剪
trim PE SRR519926_1.fastq good1.fq bad1.fq SRR519926_2.fastq  good2.fq bad2.fq SLIDINGWINDOW:4:25 MINLEN:36 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10
fastqc good1.fq
```


