## PacBio入门


### 工具准备

```bash
#更新stuff
sudo apt-get update

#基本软件安装
sudo apt-get -y install screen git curl gcc make g++ python-dev unzip \
default-jre pkg-config

#安装PBcR所需的perl模块
sudo cpan App::cpanminus
sudo cpan Statistics::Descriptive

#安装wgs-assembler
wget http://sourceforge.net/projects/wgs-assembler/files/wgs-assembler/wgs-8.2beta/wgs-8.2beta-Linux_amd64.tar.bz2
tar -jxf wgs-8.2beta-Linux_amd64.tar.bz2

#将wgs添加到$PATH
PATH=$PATH:$HOME/wgs-8.2beta/Linux-amd64/bin/
```

### 数据准备

下载样例Lambda phage数据集（很小而可以快速地进行处理）。更多挑战性的测试（大数据大电脑）见PacBio数据集：https://github.com/PacificBiosciences/DevNet/wiki/Datasets

```bash
#Download the sample data
wget http://www.cbcb.umd.edu/software/PBcR/data/sampleData.tar.gz
tar -zxf sampleData.tar.gz
cd sampleData/
```
### 数据处理

将FastA转换为faux-fastQ:

```bash
#This is really old PacBio data, provided in fastA format. 
#Look at the reads - note that they are not actually as long as 
#I just told you they should be. The PacBio tech has improved massively over the past few years.

java -jar convertFastaAndQualToFastq.jar \
pacbio.filtered_subreads.fasta > pacbio.filtered_subreads.fastq
```

使用wgs进行assembly：

```bash
PBcR -length 500 -partitions 200 -l lambda -s pacbio.spec \
-fastq pacbio.filtered_subreads.fastq genomeSize=50000
```

查看输出，可以发现phage基因组被分成2个contig。
