# 2-MeRIP-seq基本下游分析

作者 · 不加酵母的面团（DoughWithoutYeast）

原发表于 2026-10-27 · [博客园原文](https://www.cnblogs.com/DoughWithoutYeast/p/23450577)

## 1\. HOMER 寻找/富集motif

主要参考的博客 `https://www.jianshu.com/p/467d970ec097`\
博客内非常详尽地展示了HOMER的安装和不同模式的基本使用.\
对于MeRIP-seq的数据来说，我们一般使用`findMotifsGenome.pl`寻找peak区域的motif。

### 1.1 正是HOMER输入bed文件,下载HOMER对应基因组

首先我们需要构建适合HOMER的输入文件。输入文件要求5列。参考下面的格式。可以在R或者是python内将bed文件或者是peak文件进行简单整理即可！

```
chr1    6240163 6240342 T1_peak_4_peak_1     +
chr1    6240163 6240342 T1_peak_1_peak_2     +
chr1    6248764 6249158 T1_peak_6_peak_1     +
chr1    6248764 6249158 T1_peak_3_peak_2     +
chr1    6249670 6249772 T1_peak_7_peak_1     +
```

随后我们下载HOMER配套的基因组。

```bash
cd /path/to/homer
configureHomer.pl -install mm10
```

### 1.2 正式HOMER

```bash
findMotifsGenome.pl T_for_HOMER.bed -mask mm10 ./homer/T_motifDir -size 200 -len 6 -p 12 -rna
```

* `-mask mm10`：屏蔽掉mm10版本上的重复序列，也就是咱们上面下载得mm10
* `size 200`：会以每个 peak 的中心点为基准，向两侧各延伸 100 bp，从而提取一个总长 200 bp 的序列窗口进行 motif 分析。
* `-len`：指定要寻找的 motif 长度
* `-p`：运行使用的核数
* `-rna`：HOMER 会将序列视为 RNA，在 motif 搜索中使用尿嘧啶（U）。m6A是RNA上的修饰，如果不加`-rna`，HOMER会用T去匹配，导致无法正确识别以`U`结尾的m6A motif（如GGACU），从而显著降低发现真实motif的灵敏度。

## 2\. diffbind进行差异peak分析

diffbind是专用于分析差异peak的R包。参考`https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf` \& `https://github.com/nshanian/DiffBind`\
Diffbind主要依托edgR和DESeq2.它做的事情是，根据bam文件和peak文件，对每一个peak上的reads进行定量之后，进行差异表达分析，只是主体不是gene而是peak。输出的结果也是类似用featurecount定量后再使用edgR/DESeq2的结果。带有LFC和adjp。

### 2.1 构建config csv文件

```r
config <- data.frame(SampleID = c('C1','C2','T1','T2'),
                     Tissue = c(NA,NA,NA,NA),
                     Factor = c(NA,NA,NA,NA),
                     Condition = c('C','C','T','T'),
                     Treatment = c(NA,NA,NA,NA),
                     Replicate = c(1,2,1,2),
                     bamReads = c(
                       '/path/to/C1.m6A_sorted.markdup.bam',
                       '/path/to/C2.m6A_sorted.markdup.bam',
                       '/path/to/T1.m6A_sorted.markdup.bam',
                       '/path/to/T2.m6A_sorted.markdup.bam'
                     ),
                     ControlID = c('C1','C2','T1','T2'),
                     bamControl = c(
                       '/path/to/C1.input_sorted.markdup.bam',
                       '/path/to/C3.input_sorted.markdup.bam',
                       '/path/to/T1.input_sorted.markdup.bam',
                       '/path/to/T2.input_sorted.markdup.bam'
                     ),
                     Peaks = c(
                       '/path/to/C1_peaks.narrowPeak',
                       '/path/to/C2_peaks.narrowPeak',
                       '/path/to/T1_peaks.narrowPeak',
                       '/path/to/T2_peaks.narrowPeak'
                     ),
                     PeakCaller = c('narrow','narrow','narrow','narrow'))
config$Peaks <- as.character(config$Peaks)

write.csv(config,file = './config.csv',row.names = F)
```

### 2.2 run diffbind

```r
library(DiffBind)

dbObj   <- dba(sampleSheet="./config.csv")
dbObj <- dba.count(dbObj, bUseSummarizeOverlaps=TRUE)
# Establishing a contrast
dbObj <- dba.contrast(dbObj,minMembers = 2)
dbObj <- dba.analyze(dbObj, method=DBA_ALL_METHODS)
# summary of results
dba.show(dbObj, bContrasts=T)
# overlapping peaks identified by the two different tools (DESeq2 and edgeR)
dba.plotVenn(dbObj,contrast=1,method=DBA_ALL_METHODS)
comp1.edgeR <- dba.report(dbObj, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(dbObj, method=DBA_DESEQ2, contrast = 1, th=1)
# EdgeR
out <- as.data.frame(comp1.edgeR)
write.csv(out, file="./edgeR.csv", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.csv(out, file="./deseq2.csv", sep="\t", quote=F, col.names = NA)
# Create bed files for each keeping only significant peaks (p < 0.05)
# EdgeR
out <- as.data.frame(comp1.edgeR)
edge.bed <- out[ which(out$FDR < 0.05), 
                 c("seqnames", "start", "end", "strand", "Fold")]
edge.bed$seqnames <- paste0('chr',edge.bed$seqnames)
write.table(edge.bed, 
            file = './edgeR_sig.bed', 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
# DESeq2
out <- as.data.frame(comp1.deseq)
deseq.bed <- out[ which(out$FDR < 0.05), 
                  c("seqnames", "start", "end", "strand", "Fold")]
deseq.bed$seqnames <- paste0('chr',deseq.bed$seqnames)
write.table(deseq.bed, 
            file = './deseq2_sig.bed', 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE, 
            col.names = FALSE)
```

输出的数据长这样的：

```
     X seqnames     start       end width strand      Conc   Conc_T    Conc_C      Fold      p.value        FDR
1 5422        2 119947821 119948221   401      *  8.979577 4.102712  9.954816 -4.248766 8.091028e-06 0.01022994
2 9436        8  79673424  79673824   401      * 10.220829 7.052768 11.138257 -2.486556 1.191287e-05 0.01022994
3 8026        6  95121873  95122571   699      *  8.034132 1.013865  9.028564 -6.372807 1.571226e-05 0.01022994
4 2682       13  97126598  97126998   401      *  8.217912 2.327993  9.205696 -5.262004 1.595063e-05 0.01022994
5 5236        2  52112151  52112551   401      *  8.788550 3.098195  9.774512 -5.067777 1.697395e-05 0.01022994
6 8893        7 118134883 118135283   401      *  9.882239 6.127048 10.827803 -3.101252 1.801783e-05 0.01022994
```

是不是同差异表达基因分析的输出结果很像嘞？

其实diffbind内的教程内有很多关于如何分析ChIP-seq-liked数据的角度和教程，我觉得都是非常有意义的，这里就不赘述，感兴趣的可以移步看看！

## 3\. 使用input做bulk RNA-seq与peak located gene overlap

在上游分析时，我们提过一嘴，MeRIP-seq其实同其他ChIP-seq一样，都是对input进行高深度测序，然后再测一个IP。因此我们可以认为input其实就是同一般的bulk RNA-seq一样，一般可以用input来做RNA-seq分析，看那些gene在转录水平即表达发生了该边。如果同m6A peak进行overlap可以得到更加多维的观察视角。

| type | 表达上调 | 表达下调 | 表达不变 |
| --- | --- | --- | --- |
| m6A⬆ | 有可能m6A增加是由于该基因表达上升改变引起的 | m6A上升幅度被低估 | m6A上升 |
| m6A⬇ | m6A下降幅度被低估 | m6A下降可能是由于基因表达下降引起的 | m6A下降 |

## 4.定量m6A水平

定量的思路可以参考使用TPM定量gene表达的思路。其实peak同gene也无分别，都是基因组上的一段区域，只是peak相对较小而gene较大。同时，m6A定量，或者是其他的ChIP-seq的数据定量peak，计算的往往不是一个peak的"TPM"，而是FPKM ratio。IP/input的ratio。

### 4.1 get overlap peak

取出每组内生物学重复的overlap的peak

```bash
bedtools intersect -a /path/to/C1_peaks.narrowPeak  -b /path/to/C2_peaks.narrowPeak -nonamecheck >C.bed
bedtools intersect -a /path/to/T1_peaks.narrowPeak -b /path/to/T2_peaks.narrowPeak -nonamecheck >T.bed
```

### 4.2 get coverage

对overlap的peak merge后统计bam文件在peak上的coverage

```bash
cat T.bed C.bed | sort -k1,1 -k2,2n | mergeBed -i - | bedtools multicov \
-bams \
/path/to/T1.m6A_sorted.markdup.bam \
/path/to/T1.input_sorted.markdup.bam \
/path/to/T2.m6A_sorted.markdup.bam \
/path/to/T2.input_sorted.markdup.bam \
/path/to/C1.m6A_sorted.markdup.bam \
/path/to/C1.input_sorted.markdup.bam \
/path/to/C2.m6A_sorted.markdup.bam \
/path/to/C2.input_sorted.markdup.bam \
-bed -> Merge.coverage
```

### 4.3 get bam total reads

统计每个bam文件对应的reads，用来归一化

```bash
samtools idxstats /path/to/T1.m6A_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
12833759
samtools idxstats /path/to/T1.input_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
26811273
samtools idxstats /path/to/T2.m6A_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
7167560
samtools idxstats /path/to/T2.input_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
18028653
samtools idxstats /path/to/C1.m6A_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
22882246
samtools idxstats /path/to/C1.input_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
20611004
samtools idxstats /path/to/C2.m6A_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
28676000
samtools idxstats /path/to/C2.input_sorted.markdup.bam | cut -f3 | awk 'BEGIN {total=0}{total += $1}END {print total}'
18359125
```

### 4.4 get FPKM

用各自的total reads计算FPKM

```bash
awk 'OFS="\t" {print $0, \
($4*10^9)/(($3-$2)*12833759), \
($5*10^9)/(($3-$2)*26811273), \
($6*10^9)/(($3-$2)*7167560), \
($7*10^9)/(($3-$2)*18028653), \
($8*10^9)/(($3-$2)*22882246), \
($9*10^9)/(($3-$2)*20611004), \
($10*10^9)/(($3-$2)*25428676), \
($11*10^9)/(($3-$2)*16859648)}' Merge.coverage > Merge.FPKM
```

### 4.5 get ratio

`ratio = IP / input`

```bash
awk 'OFS="\t" {print $0,(($12+1)/($13+1)),\
(($14+1)/($15+1)),\
(($16+1)/($17+1)),\
(($18+1)/($19+1))}' Merge.FPKM > Merge.ratio
```

用上面的ratio就可以基本比较不同组别之间的m6A水平了，但是具体检验还是参考diffbind比较合适！

---

[合集目录](README.md) · [上一篇](01-upstream-analysis.md)
