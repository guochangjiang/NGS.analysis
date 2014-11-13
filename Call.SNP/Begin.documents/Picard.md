## Picard工具箱

### 1. 简介
Picard由一系列用于处理BAM格式相关的下代测序数据的工具组成。

Picard命令行工具打包成可执行的jar文件，用法如下：

	java jvm-args-jar PicardCommand.jar OPT1=value1 OPT2=value2...

大多数命令设计于2G的JVM，所以推荐使用JVM参数`-Xmx2g`。

### 2. 标准选项

大多数Picard程序都有一下选项：

|选项|描述|
|--------|---------|
|--help|展示该工具的特殊选项|
|--stdhelp|展示该工具的特殊选项和其他所有Picard命令行工具的共有选项|
|--version|展示程序版本|
|TMP_DIR=File|该选项可指定0次或更多次|
|VERBOSITY=LogLevel|控制记录的冗余度，默认值为INFO。可以设定为null来清除默认值。可能的值为ERROR, WARNING,INFO,DEBUG.|
|QUIET=Boolean|是否禁止输出工作总结信息之系统err。默认值为false.该选项可以设定为null以清除默认值。可能的值有true, false.|
|VALIDAYTION_SRTRINGENCY=ValidationStringency|对有该程序读取的所有SAM文件read的验证严谨度。将严谨度设为SILENT可以提高处理BAM文件的绩效，此时可变长度的数据，如read, qualities, tags等不需要进行解码。默认值为STRICT。可将该选项设定为null以清除默认设置。可能值为STRICT, LENIENT, SILENT.|
|COMPRESSION_LEVEL=Integer|所有创建的压缩文件的压缩程度（如BAM和GELI）.默认值为5.可将该选项设定为null以清除默认值。|
|MAX_RECORDS_IN_RAM=Integer|当写需sort的SAM文件时，该选项可以指定在写到到磁盘前保存在RAM中的记录数。增加该值可以减少文件操作次数但提高RAM占有。默认值为500,000，该选项可设定为null以清除默认设置。|
|CREATE_INDEX=Boolean|当写一个coordinate-sorted的BAM文件时是否创建BAM索引。默认值为false。该选项可设定为null以清除默认设置。可能的取值：true, false.|
|CREATE_MD5_FILE=Boolean|是否为BAM或FASTAQ文件创建MD5码|

### 3. 工具集

#### 3.1 AddCommentsToBam
向指定的BAM文件的头部增加更多内容。将修改过头部的文件拷贝到指定输出文件。成块拷贝方法的采取以保证有效的转移到输出文件。不支持SAM文件。

|选项|描述|
|----|----|
|INPUT=File|指定输入文件，**必需**|
|OUTPUT=FIle|指定输出文件，**必需**|
|COMMENT=String|指定需要添加到BAM文件头的内容，可多次指定|

#### 3.2 AddOrReplaceReadGroups
用新的read group取代输出文件中的所有read group，并将输出文件中的所有read分配到该read group

|选项|描述|
|----|----|
|INPUT=File|指定输入文件，**必需**|
|OUTPUT=FIle|指定输出文件，**必需**|
|SORT_ORDER=SortOrder|选择输出的sort order。如不指定则与输入文件相同。默认值为null。可取值：unsorted, queryname, coordinate.|
|RGID=String|read group的默认值: 该值可以设定为null以清除默认值|
|RGLB=String|read group library，**必需**.|
|RGPL=String|read group platform (如illumina, solid)，**必需**|
|RGPU=String|read group platform unit (如run barcode) ，**必需**|
|RGSM=String|read group sample name,**必需**|
|RGCN=String|read group sequencing center name,默认值：null.|
|RGDS=String|read group description，默认值: null.|
|RGDT=Iso8601Date|read group run date，默认值: null.|
|RGPI=Integer|read group predicted insert size，默认值: null.|

#### 3.3 BamToBfq
用法：`BamToBfq [options]

为Maq aligner创建BFQ文件。

|选项|描述|
|----|----|
|INPUT=File|要解析的BAM文件，**必需**|
|ANALYSIS_DIR=File|二进制输出文件的分析目录，**必需**|
|FLOWCELL_BARCODE=String|流动细胞条形码,如30PYMAAXX,**必需**，不能和选项OUTPUT_FILE_PREFIX一起使用|
|LANE=Integer|Lane数目，默认值：null。不能和选项OUTPUT_FILE_PREFIX一起使用|
|OUTPUT_FILE_PREFIX=String|所有输出文件的前缀，**必需**。不能和上两项同时使用。|
|READS_TO_ALIGN=Integer|align的read数量，默认：null（null=all)|
|READ_CHUNK_SIZE=Integer|用于分入单独group的read数目，默认值2,000,000。该选项可以设置为null以清除默认设置。|
|PAIRED_RUN=Boolean|是否是pair-end，**必须**。可能取值：true, false|
|RUN_BARCODE=String|弃用选项，用READ_NAME_PREFIX而非默认值null，不可和选项READ_NAME_PREFIX同时使用|
|READ_NAME_PREFIX=String|从每个read名中剥去的前缀，以使满足Maq的要求。默认值：null。|
|INCLUDE_NON_PF_READS=Boolean|是否保留non-PF read，默认值为false.该选项可设定为null以清除默认设置。可能取值：true, false.|
|CLIP_ADAPTERS=Boolean|是否从read剪掉adapter，默认值为true。可设定为null以清除默认设置。可能取值：true, false.|
|BASES_TO_WRITE=Integer|从每个read的那些碱基写到BFQ文件。如为non-null，只有每个read的前BASES_TO_WRITE个碱基会写到bfq文件。默认值为null|

#### 3.4 BamIndexStats

