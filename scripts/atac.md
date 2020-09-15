
## Contents

- [Libraries](atac.md#libraries)
- [Software and operating system](atac.md#software-and-operating-system)
- [Data processing](atac.md#data-processing)
  - [Quality check](atac.md#quality-check)
  - [Trim adapters and filter by base quality](atac.md#trim-adapters-and-filter-by-base-quality)
- [Alignment](atac.md#alignment)
  - [Prepare and index references](atac.md#prepare-and-index-references)
  - [Align and sort](atac.md#align-and-sort)
  - [Mark duplicates, index and flagstat](atac.md#mark-duplicates-index-and-flagstat)
  - [Filter and index](atac.md#filter-and-index)
  - [Check chrM abundance](atac.md#check-chrM-abundance)
  - [Check number of human and virus reads, and depth](atac.md#check-number-of-human-and-virus-reads-and-depth)
  - [Insert size plots](atac.md#insert-size-plots)
  - [Create bigwig files](atac.md#create-bigwig-files)
  - [ATACseqQC](atac.md#atacseqqc)
- [Differential analysis](atac.md#differential-analysis)
  - [Define mappable genome](atac.md#define-mappable-genome)
  - [Obtain human mappable regions](atac.md#obtain-human-mappable-regions)
  - [Coverage in regions](atac.md#coverage-in-regions)
  - [Testing](atac.md#testing)



## Libraries

id | files
---|------
WTnoVi | WTnoVi.R1.fastq.gz,WTnoVi.R2.fastq.gz
WTnoVLP | WTnoVLP.R1.fastq.gz,WTnoVLP.R2.fastq.gz
WTConVLP | WTConVLP.R1.fastq.gz,WTConVLP.R2.fastq.gz
WTVprVLP | WTVprVLP.R1.fastq.gz,WTVprVLP.R2.fastq.gz

Saved in folder `~/fastq_atac`



## Software and operating system

Software:

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Standard Unix tools: awk, cut, grep, pigz, sort, echo, mkdir
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [bwa v0.7.15-r1140](http://bio-bwa.sourceforge.net/)
- [samtools v1.3.1](http://samtools.sourceforge.net/)
- [sambamba v0.6.5](https://academic.oup.com/bioinformatics/article/31/12/2032/214758)
- [bedtools v2.27.0](http://bedtools.readthedocs.io/en/latest/)
- [deeptools v2.4.2-5-f439d22](https://deeptools.readthedocs.io/en/develop/index.html)
- [ATACseqQC v1.6.4](https://www.bioconductor.org/packages/release/bioc/html/ATACseqQC.html)
- [R v3.3.2](https://www.r-project.org/). Libraries:
  - [ggplot2 v2.2.1](http://ggplot2.org/)
  - [data.table v1.10.4](https://cran.r-project.org/web/packages/data.table/index.html)

Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)



## Data processing

### Quality check

```bash
cd ~/fastq_atac

mkdir ../fastqc_atac

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc_atac/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc_atac $fq"
done
```


### Trim adapters and filter by base quality

```bash
cd ~/fastq_atac

mkdir ../fastq_trimmed_atac

for fq1 in *.R1.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%.R1.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed_atac/$bname.log --mem 4G --wrap "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 -o ../fastq_trimmed_atac/$fq1 -p ../fastq_trimmed_atac/$fq2 $fq1 $fq2 > ../fastq_trimmed_atac/$bname.txt"
done
```



## Alignment

### Prepare and index references

Downloaded human hg38 reference genome from [iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) and HIV reference genome `HIV1_1LTR_circle.fa` from the [genomes folder](../genomes/HIV1_1LTR_circle.fa)

```bash
cd ~/reference
cat hg38.fa HIV1_1LTR_circle.fa > hg38_HIV1_1LTR_circle.fa

bwa index hg38_HIV1_1LTR_circle.fa
samtools faidx hg38_HIV1_1LTR_circle.fa

rm hg38.fa
rm HIV1_1LTR_circle.fa
```


### Align and sort

```bash
cd ~/fastq_trimmed_atac

mkdir ../bam_atac

ref=~/reference/hg38_HIV1_1LTR_circle.fa

for fq1 in *.R1.fastq.gz
do
  bname=${fq1%.R1.fastq.gz}
  fq2=${fq1/_R1_/_R2_}
  sbatch -J $bname -o ../bam_atac/$bname.log --mem 32G --wrap "bwa mem -t 20 -M $ref $fq1 $fq2 | \
  samtools view -@ 20 -b - | \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o ../bam_atac/$bname.tmp.bam -"
done
```


### Mark duplicates, index and flagstat

```bash
cd ~/bam_atac

mkdir ../flagstat_atac

for bam in *.tmp.bam
do
  bname=${bam%.tmp.bam}
  sbatch -J $bname -o $bname.markdup.log --mem 16G --wrap "sambamba markdup -t 20 $bam $bname.bam 2> $bname.markdup.txt && \
  samtools flagstat $bname.bam > ../flagstat_atac/$bname.txt"
done
```


### Filter and index

Create whitelist from blacklist:

```bash
cd ~/reference
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
pigz -d hg38-blacklist.v2.bed.gz
cut -f1-2 hg38_HIV1_1LTR_circle.fa.fai | awk -v OFS="\t" '{print $1, 0, $2}' | bedtools subtract -a - -b hg38-blacklist.v2.bed | sort -k1,1 -k2,2n > hg38-whitelist.v2.bed
grep -v "chrM" hg38-whitelist.v2.bed > hg38-whitelist.v2.nochrM.bed
```

Removing duplicate and supplementary alignments, reads mapping to chrM and blacklisted regions as well as filtering by mapping quality:

```bash
cd ~/bam_atac

for bam in `ls *.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.clean.log --mem 16G --wrap "samtools view -b -L ~/reference/hg38-whitelist.v2.nochrM.bed -q 10 -F 3844 -@ 20 $bam > $bname.clean.bam && \
  samtools index $bname.clean.bam"
done
```

Check number of aligned reads:

```bash
cd ~/bam_atac

for bam in `ls *.bam | grep -v "tmp"`
do
  n=`samtools view -@ 20 $bam -c`
  echo $bam, $n
done
```


### Check chrM abundance

```bash
cd ~/bam_atac

for bam in `ls *.bam | grep -v "clean"`
do
  n_tot=`samtools view -@ 20 $bam -c`
  n_chrM=`samtools view -@ 20 $bam chrM -c`
  pct=`echo "scale=1; 100*$n_chrM/$n_tot" | bc`
  echo -e "$bam\t$n_tot\t$n_chrM\t$pct%"
done
```


### Check number of human and virus reads, and depth

```bash
cd ~/bam_atac

for bam in *.clean.bam
do
  n_human=`samtools idxstats $bam | grep chr* | cut -f1 | xargs samtools view -@ 20 $bam -c`
  depth_human=`samtools depth -r chr1 $bam | awk '{sum+=$3} END {print sum/NR}'`
  n_virus=`samtools view -@ 20 $bam HIV1_1LTR_circle -c`
  depth_virus=`samtools depth -r HIV1_1LTR_circle $bam | awk '{sum+=$3} END {print sum/NR}'`
  echo -e "$bam\t$n_human\t$depth_human\t$n_virus\t$depth_virus"
done
```


### Insert size plots

```bash
cd ~/bam_atac

mkdir ../stats_atac

for bam in *.clean.bam
do
  bname=${bam%.bam}
  mkdir ../stats_atac/$bname
  sbatch -J $bname -o ../stats_atac/$bname/$bname.log --mem 8G --wrap "samtools stats $bam > ../stats_atac/$bname/$bname.txt && \
  ~/sw/samtools/samtools-1.3.1/bin/plot-bamstats -p ../stats_atac/$bname/ ../stats_atac/$bname/$bname.txt"
done
```


### Create bigwig files

Standard approach:

```bash
cd ~/bam_atac

mkdir ../bw_atac

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bw_atac/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done
```

With qPCR fold-change normalisation using the `--scaleFactor` option::

```bash
cd ~/bam_atac

for bam in *.clean.bam
do
  bname=${bam%.bam}
  echo $bam
  samtools view -b -@ 20 $bam HIV1_1LTR_circle > ${bname}.HIV1_1LTR_circle.bam && samtools index ${bname}.HIV1_1LTR_circle.bam  
done

# WTnoVi
bam=WTnoVi.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1.17 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# WTnoVLP
bam=WTnoVLP.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# WTVprVLP
bam=WTVprVLP.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# WTConVLP 
bam=WTConVLP.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1.5 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &
```


### ATACseqQC

```r
library(ATACseqQC)

# estimate library complexity
f <- "WTnoVi.clean.bam"
freq <- readsDupFreq(f)
estimateLibComplexity(freq)

# fragment size distribution
#dir.create("../atacseqqc/frag_size_dist")

for (f in dir(pattern = "*clean.bam$")){
  bname <- unlist(strsplit(f, "_"))[1]
  print(bname)
  pdf(paste("../atacseqqc/frag_size_dist/", bname, ".pdf", sep = ""))
  fragSizeDist(f, bname)
  dev.off()
}
```



## Differential analysis

### Define mappable genome

```bash
cd ~/bam_atac

g=~/reference/hg38_HIV1_1LTR_circle.fa.fai

sbatch -J mappable -o mappable.log --mem 16G --wrap "bedtools genomecov -ibam WTnovi_S1.clean.bam -bg | sort -k1,1 -k2,2n | mergeBed > WTnoVi.clean.bed"
```


### Obtain human mappable regions

Generate 100,000 mappable regions at random from human each the size of the viral genome, then add viral genome:

```bash
cd ~/bam_atac

mkdir ../de_atac

g=~/reference/hg38_HIV1_1LTR_circle.fa.fai

bedtools intersect \
-a <(bedtools random -l 3825 -n 150000 -seed 123 -g <(cut -f1,2 $g | grep "^chr" | grep -E -v "random|chrUn|chrEBV|chrM")) \
-b WTnoVi.clean.bed \
-wa -u | cut -f1-3 | head -n 100000 | bedtools sort -i > ../de_atac/regions.bed

cd ../de_atac
echo -e 'HIV1_1LTR_circle\t1\t3826' >> regions.bed
```


### Coverage in regions

```bash
cd ~/bam_atac

sbatch -J multicov -o ../de_atac/multicov.log --mem 16G --wrap "bedtools multicov -bams WTnoVi.clean.bam WTnoVLP.clean.bam WTConVLP.clean.bam WTVprVLP.clean.bam -bed ../de_atac/regions.bed > ../de_atac/regions_counts.bed"
```


### Testing

```r
#cd ~/de_atac
#R

library(data.table)
library(ggplot2)

# Enlarge the view width when printing tables
options(width = 250)

# Load data
data <- fread("regions_counts.bed")
setnames(data, c("chr", "start", "end", "WTnoVi", "WTnoVLP", "WTConVLP", "WTVprVLP"))


########################
# WTVprVLP vs WTConVLP #
########################
detable_VprVLPvi_ConVLPvi <- copy(data[, c("chr", "start", "end", "WTConVLP", "WTVprVLP")])

detable_VprVLPvi_ConVLPvi[chr == "HIV1_1LTR_circle", WTVprVLP := round(WTVprVLP/30)]
detable_VprVLPvi_ConVLPvi[chr == "HIV1_1LTR_circle", WTConVLP := round(WTConVLP/30)]

n_WTVprVLP <- sum(detable_VprVLPvi_ConVLPvi$WTVprVLP)
n_WTConVLP <- sum(detable_VprVLPvi_ConVLPvi$WTConVLP)

detable_VprVLPvi_ConVLPvi[, log2fc := log2((WTVprVLP/n_WTVprVLP)/(WTConVLP/n_WTConVLP))]
detable_VprVLPvi_ConVLPvi[, pval := fisher.test(matrix(c(WTVprVLP, WTConVLP, n_WTVprVLP, n_WTConVLP), nrow = 2), alternative = "two.sided")$p.value, by = 1:nrow(detable_VprVLPvi_ConVLPvi)]
detable_VprVLPvi_ConVLPvi[, pvaladj := p.adjust(detable_VprVLPvi_ConVLPvi$pval, method = "BH")]

# Volcano plot
gg <- ggplot(detable_VprVLPvi_ConVLPvi[is.finite(log2fc)], aes(x = log2fc, y = -log10(pvaladj))) +
geom_point(size = 1, alpha = 0.25) +
geom_point(data = detable_VprVLPvi_ConVLPvi[chr == "HIV1_1LTR_circle"], aes(x = log2fc, y = -log10(pvaladj)), color='red', size=3) +
geom_point(data = detable_VprVLPvi_ConVLPvi[chr == "chr13" & start == 60663053 & end == 60666878], aes(x = log2fc, y = -log10(pvaladj)), color='green', size=3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
theme_classic() +
coord_cartesian(xlim = c(-4, 4), clip = "off") +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"))

ggsave('../figures/VprVLPvi_ConVLPvi_volcano.png')


# Tables
write.table(detable_VprVLPvi_ConVLPvi[order(pvaladj)], file = "~/tables/VprVLPvi_ConVLPvi.30xnorm.only.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)


#######################
# ConVLPvi vs WTnoVLP #
#######################
detable_ConVLPvi_WTviB <- copy(data[, c("chr", "start", "end", "WTnoVLP", "WTConVLP")])

detable_ConVLPvi_WTviB[chr == "HIV1_1LTR_circle", WTConVLP := round(WTConVLP/30)]
detable_ConVLPvi_WTviB[chr == "HIV1_1LTR_circle", WTnoVLP := round(WTnoVLP/30)]

n_WTConVLP <- sum(detable_ConVLPvi_WTviB$WTConVLP)
n_WTnoVLP <- sum(detable_ConVLPvi_WTviB$WTnoVLP)

detable_ConVLPvi_WTviB[, log2fc := log2((WTConVLP/n_WTConVLP)/(WTnoVLP/n_WTnoVLP))]
detable_ConVLPvi_WTviB[, pval := fisher.test(matrix(c(WTConVLP, WTnoVLP, n_WTConVLP, n_WTnoVLP), nrow = 2), alternative = "two.sided")$p.value, by = 1:nrow(detable_ConVLPvi_WTviB)]
detable_ConVLPvi_WTviB[, pvaladj := p.adjust(detable_ConVLPvi_WTviB$pval, method = "BH")]

# Volcano plot
gg <- ggplot(detable_ConVLPvi_WTviB[is.finite(log2fc)], aes(x = log2fc, y = -log10(pvaladj))) +
geom_point(size = 1, alpha = 0.25) +
geom_point(data = detable_ConVLPvi_WTviB[chr == "HIV1_1LTR_circle"], aes(x = log2fc, y = -log10(pvaladj)), color='red', size=3) +
geom_point(data = detable_ConVLPvi_WTviB[chr == "chr2" & start == 22400721 & end == 22404546], aes(x = log2fc, y = -log10(pvaladj)), color='green', size=3) +
ylab(expression("-log"[10]*"FDR")) +
xlab(expression("log"[2]*"FC")) +
theme_classic() +
coord_cartesian(xlim = c(-4, 4), clip = "off") +
theme(axis.title = element_text(size=14), axis.text.y = element_text(size=14, color = "black"), axis.text.x = element_text(size=14, , color = "black"))

ggsave('~/figures/ConVLPvi_WTviB_volcano.png')


# Tables
write.table(detable_ConVLPvi_WTviB[order(pvaladj)], file = "~/tables/ConVLPvi_WTviB.30xnorm.only.txt", row.names = FALSE, col.names = TRUE, sep = '\t', quote = FALSE)
```
