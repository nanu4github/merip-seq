# 1-MeRIP-seq基础上游分析

作者 · 不加酵母的面团（DoughWithoutYeast）

原发表于 2026-10-26 · [博客园原文](https://www.cnblogs.com/DoughWithoutYeast/p/23442784)

## 0\. overall

我们先来整体看看当我们拿到一个MeRIP-seq的数据，不管是自己测得还是从网络上下载的，我们需要进行什么分析呢\
主要包括：

1. 质量控制：文件完整性\+测序质量评估
2. 去接头，顺带去测序质量不佳的reads
3. 去rRNA
4. 比对到参考基因组
5. 去PCR重复
6. call peak

以上是我们的上游分析的结果。下游可能还有关于motif的富集分析，使用diffbind进行差异peak的分析。敬请关注一下篇！

## 1\. 质控

### 1.1 检查文件完整性

**每次都要强调，不管是自己的测序数据还是从网络上下载的数据，第一步都是检查文件完整性！！！**\
`vim check_md5.sh`将下面的内容粘贴至该脚本内，保存后直接运行`bash check_md5.sh`。\
该文件会遍历当前目录下面的fastq gz文件，并读取对应的同门的md5文件，检查gz文件的md5码是否于md5文件内的内容一致。

```bash
#!/bin/bash
# 遍历当前目录下所有 .fq.gz 文件
for fq_file in *.fq.gz; do
    # 检查对应的 .md5 文件是否存在
    md5_file="${fq_file}.md5"
    if [[ -f "$md5_file" ]]; then
        # 计算当前 .fq.gz 文件的 MD5 值
        computed_md5=$(md5sum "$fq_file" | awk '{print $1}')
        # 读取 .md5 文件的内容
        file_md5=$(head -n 1 "$md5_file" | awk '{print $1}')
        # 比较两者是否一致
        if [[ "$computed_md5" == "$file_md5" ]]; then
            echo "√√√√√√√√√√√√√√√√√√√√√√√√ MD5 matches for $fq_file"
        else
            echo "×××××××××××××××××××××××× MD5 mismatch for $fq_file"
            echo "Computed: $computed_md5"
            echo "Expected: $file_md5"
        fi
    else
        echo "！！！！！！！！！！！！！！ MD5 file missing for $fq_file"
    fi
done
```

### 1.2 fastqc \+ multiqc 检查测序数据质量

```bash
#!/bin/bash
cd /path/to/data_path/ #进入到fastq.gz文件保存目录下
ls ./*fq.gz | xargs fastqc -t 12 -o  /path/to/01qc/  #输出在01qc目录下       
source ~/.bashrc #先激活以下环境以便使用各种软件等
conda activate base
multiqc /path/to/01qc/ 
```

### 1.3 解读质控结果

* 首先需要满足所有的数据文件完整之后才可以进行包括质控在内的所有分析
* 看multiqc的整体多样品报告时，可能会出现一下报错，例如GC content的报错，这是MeRIP-seq测序的内在特性。因为IP时拉下来的片段可能就是具有某种特征的序列，因此，不仅仅是MeRIP-seq，任何ChIP-seq都有可能存在这样的特性。我们主要关注Q30等测序错误相关的指标

## 2\. 去接头

我们这里主要使用`trim_galore`去接头。

```bash
#!/bin/bash
date
trimmed_path='/path/to/02cut_adapt/trimmed_output/' #变量赋值语句
file_path='/path/to/data_path/'
for i in C1.m6A C1.input C2.m6A C2.input T1.input T1.m6A T2.input T2.m6A  #这里的每一个i都可以对应匹配到各自的fq.gz文件，用空格分开
do
        echo -e "\n ################################# \n"
        echo $file_path/${i}_1.fq.gz
        echo -e "\n ################################# \n"
        trim_galore -j 16 -q 25 --phred33 --length 36 --stringency 3 --path_to_cutadapt /path/to/cutadapt --paired $file_path/${i}_1.fq.gz $file_path/${i}_2.fq.gz --gzip -o $trimmed_path
done
echo -e  " \n ###################  DONE  ####################!!! \n "
date
```

* `-q`：设置测序质量的阈值，当测序质量低于该阈值时，将该reads丢弃。我习惯用25，当严格时可以用30
* `-j`：设置并行时使用的CPU的核数.需要根据CPU的具体情况进行改变

## 3\. 去rRNA

很多实验在提RNA时就会使用去rRNA的试剂盒对提取到的RNA进行纯化。但是我习惯再去除一遍(数据角度)\
`https://zhuanlan.zhihu.com/p/329967121` 主要参考的blog，下面使用的rRNA的index的建立等细节可以参考该blog！

```bash
#!/bin/bash
date
input_path='/path.to/02cut_adapt/trimmed_output' #变量赋值语句
output_path='/path/to/03remove_rRNA/result'
for i in C1.m6A C1.input C2.m6A C2.input T1.input T1.m6A T2.input T2.m6A
do
        echo -e "\n ################################# \n"
        echo $input_path/${i}_1_val_1.fq.gz
        echo -e "\n ################################# \n"
        bowtie2 --very-sensitive-local --no-unal -I 1 -X 1000 -p 24 -x /path/to/mm38/remove_rRNA_index/rRNA \
                -1 $input_path/${i}_1_val_1.fq.gz \
                -2 $input_path/${i}_2_val_2.fq.gz \
                -S $output_path/${i}_Map2rRNA.sam \
                --un-conc-gz $output_path/${i}_rRNAremoved.fq.gz 2>$output_path/${i}_Map2rRNAStat.xlsx
done
echo -e  " \n ###################  DONE  ####################!!! \n "
date
```

* `--very-sensitive-local`：预设灵敏度 \+ 局部比对模式
* `--no-unal`:不把未比对上的 reads 写入 SAM 文件
* `-I` `-X`：最小插入长度，最大插入长度
* `-p`：并行是使用的核数
* `-x`：rRNA的bowtie2的index，在上面提供的blog连接内有详细介绍
* `-1` `-2`：分别输入fastq 1 \& fastq 2文件路径
* `-S`：比对到rRNA上的reads保存的sam文件
* `--un-conc-gz`：把未比对上的 reads 对（unconcordant，即没有比对到 rRNA 的配对）以 gzip 压缩写入指定文件

我们可以在定向输出的xlsx内查看比对到rRNA上的比率有多少。\
类似：

```
58963895 reads; of these:
  58963895 (100.00%) were paired; of these:
    53258887 (90.32%) aligned concordantly 0 times
    5413320 (9.18%) aligned concordantly exactly 1 time
    291688 (0.49%) aligned concordantly >1 times
    ----
    53258887 pairs aligned concordantly 0 times; of these:
      46037 (0.09%) aligned discordantly 1 time
    ----
    53212850 pairs aligned 0 times concordantly or discordantly; of these:
      106425700 mates make up the pairs; of these:
        104726907 (98.40%) aligned 0 times
        1695014 (1.59%) aligned exactly 1 time
        3779 (0.00%) aligned >1 times
11.19% overall alignment rate
```

也就是说还有11.19%是比对到rRNA上的reads，可以推断，接近10%的rRNA残留率。

## 4\. 比对

回顾一下上面的处理后我们可以拿到什么样的结果：文献完整，测序质量合格，去接头，去测序不佳的reads，去rRNA的数据(fastq format)\
我们将上面处理好的数据就可以比对到参考基因组上啦\~

```bash
#!/bin/bash
echo -e  " \n \n \n ###############  Align !!! hisat2 !!!################\n \n \n "
date
index='/path/to/hisat2_index/grcm38_snp/genome_snp'  #这里一定一定要注意物种和版本信息
pwd
input_path='path/to/03remove_rRNA/result/' 
for i in C1.m6A C1.input C2.m6A C2.input T1.input T1.m6A T2.input T2.m6A
do
  echo -e "\n ################################# \n"
  echo -e  " ${i} start mapping   "
  echo -e "\n ################################# \n"
        hisat2 -t -p 12 -x  $index \
          -1 $input_path/${i}_rRNAremoved.fq.1.gz  \
          -2 $input_path/${i}_rRNAremoved.fq.2.gz \
    -S  ./SAM/${i}.sam
  samtools sort -@ 12 -o  ./BAM/${i}_sorted.bam  ./SAM/${i}.sam
  rm ./SAM/${i}.sam
  echo -e "\n ################################# \n"
  echo -e  " ${i} done   "
  echo -e "\n ################################# \n"
done
echo -e  " \n ###################  DONE  ####################!!! \n "
date
```

hisat2部分参数

* `-t`：打印比对耗时
* `-p`：设置比对使用核数
* `x`：比对使用的index
* `-1` `-2`：比对输入的双端fastq文件
* `-S`：输出的sam文件

samtools部分参数

* `-@`：使用核数
* `-o`：输出文件

上面的代码，将比对到参考基因组的sam文件顺带收给sort了.同时还删除了较大的sam文件，保留了二进制压缩后的bam文件及其index文件。

## 5\. 去PCR重复

关于是否需要去PCR重复各家都有各自的说法。如果你觉得不必要去PCR重复的话，可以跳过这个步骤之后将后面的call peak的bam地址进行一定的修改即可。

```bash
#!/bin/bash
echo -e "\n\n\n# PCR remove duplicate !!! picard !!!\n\n\n"
date
BAM='/path/to/04mapping/BAM/'
OUT='/path/to//05remove_duplication/result/'
for i in C1.m6A C1.input C2.m6A C2.input T1.input T1.m6A T2.input T2.m6A
do
  echo -e "\n ################################# \n"
  echo ${BAM}${i}.sorted.reheader.bam
  echo -e "\n ################################# \n"
  picard MarkDuplicates \
    REMOVE_DUPLICATES=true \
    I=${BAM}${i}_sorted.bam \
    O=${OUT}${i}_sorted.markdup.bam \
    M=${OUT}${i}_markdup_metrics.txt
  echo -e "\n ################################# \n"
done
echo -e  " \n ###################  DONE  ####################!!! \n "
date
```

* `-I`：输入的bam文件
* `-O`：输出的，去重后的bam文件
* `-M`：输出的重复报告的文件

重复报告文件内含信息：

| 指标 | 含义 |
| --- | --- |
| UNPAIRED\_READS\_EXAMINED | 单端 reads 数 |
| READ\_PAIRS\_EXAMINED | 双端 read pairs 数 |
| UNMAPPED\_READS | 未比对 reads 数 |
| DUPLICATE\_PAIRS | 重复 pair 数 |
| PERCENT\_DUPLICATION | 重复率 |
| ESTIMATED\_LIBRARY\_SIZE | 估计文库复杂度 |

## 6\. call peak

最后一步啦\~\
在call peak之前的内容基本上和一般的RNA-seq无差别。\
call peak其实就是将input的样品作为背景板，看IP下来的reads，会在基因组上的哪些位置上形成peak。具体背景我将再出一期进行介绍，敬请关注！\
MeRIP-seq同其他ChIP-seq一般，在call peak时主要是使用MACS2

```bash
macs2 callpeak \
    -t /path/to/m6A/bam \
    -c /path/to/input/bam \
    -n output_name \
    -g mm \
    -B \
    --keep-dup all \
    --nomodel \
    --extsize 100 \
    --broad \
    --outdir /path/to/output_dir/
```

* `-t`：Treatment（实验组）。这里指定 m6A IP样本的BAM文件。文件必须是按坐标排序且已建立索引（.bai） 的基因组比对结果。
* `-c`：Control（对照组）。指定 Input样本的BAM文件。Input用于构建背景模型，校正RNA表达丰度和测序偏差，从而分离出真正的m6A富集信号。
* `-n`：输出文件的前缀名。MACS2会基于此生成 .narrowPeak 等结果文件。
* `-g`：有效基因组大小（Effective genome size）。
  + `mm`:是MACS2内置的简写，代表小鼠（Mouse）的有效基因组大小，数值约为 1.87e9。
  + `hs`:人类的有效基因组大小，数值约 2.7e9。
* `-B`：输出 bedGraph 文件。生成 .bdg 信号文件，便于在IGV等基因组浏览器中可视化信号强度。
* `--keep-dup`：MACS2默认会去除重复reads（--keep-dup 1）。但在MeRIP-seq中，高表达转录本上密集的reads是真实的生物学信号，默认去重会严重破坏m6A富集信号。务必加上此参数保留所有reads。尤其是在我们前面已经去除了PCR的情况下。
* `--nomodel`：RNA的片段化是随机的，不遵循DNA测序中基于链方向的位移模型。使用 --nomodel 跳过模型构建，避免得出无意义的结果。
* `--extsize`：配合`--nomodel`，手动指定延伸长度。常用的值在 100-150 bp 左右，可根据你的片段化大小调整
* `--broad`：m6A修饰在转录本上常呈现宽泛的富集区域（而非TFs结合那样的尖锐峰）
* `--outdir`：输出路径

输出的文件：

| 文件 | 格式 | 核心信息 | 主要用途 |
| --- | --- | --- | --- |
| output\_name \_peaks.xls | 文本表格 | 所有peak的汇总信息与参数 | 快速浏览、整体统计 |
| output\_name \_peaks.broadPeak | BED6\+3 | peak坐标、富集倍数、p/q值 | 下游分析的核心输入（注释、关联基因） |
| output\_name \_peaks.gappedPeak | BED12\+3 | 宽域\+窄峰的结构化信息 | 精细分析宽峰内部结构 |
| output\_name \_summits.bed | BED | peak summit坐标 | motif分析、精确位点注释 |
| output\_name \_treat\_pileup.bdg | bedGraph | IP样本的原始信号堆叠 | IGV可视化、信号检查 |
| output\_name \_control\_lambda.bdg | bedGraph | Input样本估计的局部背景 | 评估背景校正效果 |

## to be continued

以上，我们就完成了MeRIP-seq的上有分析了，下游分析常常包括motif的富集分析以及差异peak分析。这些内容将在后面的blog中补充噢\~

---

[合集目录](README.md) · [下一篇](02-downstream-analysis.md)
