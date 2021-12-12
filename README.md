# survey
Eukaryotic genome size assessment process based on second-generation sequencing


Third-party
-----------

survey package includes some third-party software:
* [python](https://www.python.org/)
* * Three-party python package
  * [pysam](https://pypi.org/project/pysam/)
  * [matplotlib](https://matplotlib.org/)
  * [numpy](https://numpy.org/doc/stable/index.html)
* [R](https://www.r-project.org/)
* [fastp](https://github.com/OpenGene/fastp)
* [fastqc](https://github.com/s-andrews/FastQC)
* [kmc](https://github.com/refresh-bio/KMC)
* [seqkit](https://github.com/shenwei356/seqkit)
* [splitfp](https://github.com/zxgsy520/splitfp)
* [SOAPdenovo](https://github.com/aquaskyline/SOAPdenovo2)
* [megahit](https://github.com/voutcn/megahit)
* [minimap2](https://github.com/lh3/minimap2)
* [samtools](https://https://github.com/samtools/samtools)
* [sam2ngs](https://github.com/zxgsy520/sam2ngs)
* [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [ghostscriptsharp](https://github.com/mephraim/ghostscriptsharp)
* [genomescope](https://github.com/schatzlab/genomescope)
* [findGSE](https://github.com/schneebergerlab/findGSE)

## Installation
```
git clone https://github.com/zxgsy520/survey.git
cd  survey/bin
chmod 755 *
cd ..
chmod 755 survey.py
cd survey
#Modify software configuration files, or add related software to environment variables.
vi config.py
```
## Quick usage
```
survey.py -h
usage: survey.py [-h] {all,ngs_qc,filter_cont,kmer_denovo} ...

URL: https://github.com/zxgsy520/survey
Survey: Eukaryotic survey analysis Pipeline

version: v2.1.0
contact: Xingguo Zhang <invicoun@foxmail.com>        

optional arguments:
  -h, --help            show this help message and exit

command:
  {all,ngs_qc,filter_cont,kmer_denovo}
    all                 all steps
    ngs_qc              quality control
    filter_cont         Sequence contamination filtering
    kmer_denovo         Genomic k-mer estimation and assembly
```
Run the entire process
```
 ./survey.py all -h
usage: survey.py all [-h] [-kl INT] [-gt {small,big,oversized}] [-d INT]
                     [-ml INT] [-m {fast,general,strict}] [-c INT] [-w INT]
                     [-s {split,}] [-mq INT] [--asm {megahit,soapdenovo,}]
                     [-q STR] -r1 FILE [FILE ...] -r2 FILE [FILE ...] -n FILE
                     [--trim INT] [--qvalue INT]
                     [--kingdom {plant,animal,fungi}] [--thread THREAD]
                     [--concurrent INT] [--refresh INT]
                     [--job_type {sge,local}] [--work_dir DIR] [--out_dir DIR]

optional arguments:
  -h, --help            show this help message and exit
  -kl INT, --kmer_length INT  #计算kmer的长度，默认17；对与一些高重复的物种，如果预估偏小可以增大该值。
                        Set the length of kmer.
  -gt {small,big,oversized}, --genome_type {small,big,oversized}   #选择基因的大小，1G以下可以使用small；1-4G可以使用big；5G以上可以使用oversized
                        Choose genome type, default=small.
  -d INT, --depth INT   Set the selected kmer depth(24, 40), default=35. #选择保留做拟合图的数据深度，一般要求24-40X（会自动选择，有时候选不准需要人工选）
  -ml INT, --maximal INT     #输出最大的kmer出现次数，一般真菌可以设置在2000（可以节约计算时间和资源）；动植物的如果是没有经验的检验不动。
                        maximal value of a counter (default=1e6).
  -m {fast,general,strict}, --mode {fast,general,strict}  #过滤污染的模型，fast为不去污染、线粒体和叶绿体；general为会去除线粒体叶绿体（污染严重也会去污染）；strict为去除污染、线粒体和叶绿体。一般的物种使用fast即可。
                        Choose the mode of operation, default=fast.
  -c INT, --cratio INT  Input pollution rate, valid when mode = strict,  #设置污染物种的占比，设置当--mode设置为general和strict有效。过滤污染占比。
                        default=10.
  -w INT, --window INT  Set the window size when calculating the GC depth.  #组装绘制GC深度图的窗口大小。
  -s {split,}, --split {split,}       #对数据进行分割并行
                        Cache output, default=       
  -mq INT, --minmapq INT      #过滤污染、线粒体和叶绿体的的比对得分。
                        Set the minimum quality value for reads
                        alignment(<20), default=0
  --asm {megahit,soapdenovo,}  #选择组装使用的软件，默认不组装。
                        Select the assembled software, if not assembled, it
                        will be empty, default=
  -q STR, --queue STR   Queue to be used for assembly tasks,default=fat.q .   #设置组装使用的队列
  -r1 FILE [FILE ...], --read1 FILE [FILE ...]      #输入二代测序的数据，支持通配符，不用对数据进行合并。
                        Second generation sequencing of reads r1.
  -r2 FILE [FILE ...], --read2 FILE [FILE ...]
                        Second generation sequencing of reads r1.
  -n FILE, --name FILE  Sample name.     #样本的名称。
  --trim INT            Set trim length, default=5   #数据修剪参数。
  --qvalue INT          The quality value that a base is qualified, default=20   #数据过滤的质量值。
  --kingdom {plant,animal,fungi}     #选择物种使用的模型
                        Choose kingdom.
  --thread THREAD       threads used in genome mapping    #设置线程

Workflow arguments:
  --concurrent INT      Maximum number of jobs concurrent (default: 10).   #设置并行的任务数目
  --refresh INT         Refresh time of log in seconds (default: 30).   #设置检测任务状况的刷新时间（30s）
  --job_type {sge,local}      #sge为调用集群运行；local为本地节点运行
                        Jobs run on [sge, local] (default: local).
  --work_dir DIR        Work directory (default: current directory).  #工作目录，分析完成后可以直接删除
  --out_dir DIR         Output directory (default: current directory). #输出结果目录
```
