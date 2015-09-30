


### 第二天：命令行运行BLAST

+ 获取BLAST：  
	curl -O ftp://ftp.ncbi.nih.gov/blast/executables/release/2.2.24/blast-2.2.24-x64-linux.tar.gz
+ 建立蛋白库  
	formatdb -i xx.protein.faa -o T -p T
+ 有筛选地BLAST  
	blastall -i query.fa -d xx.protein.faa -p blastp -b 2 -v 2 -e 1e-6 -o results.txt
+ 把BLAST输出转换为CSV格式  
   + 工具：ngs-scripts/blast-to-csv.py
   + 用法： python blast-to-csv.py results.txt > results.csv