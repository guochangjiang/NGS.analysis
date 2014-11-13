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
