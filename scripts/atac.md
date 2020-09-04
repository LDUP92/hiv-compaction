
## Contents

- [Libraries](atac.md#libraries)
- [Software and operating system](atac.md#software-and-operating-system)
- [Data processing](atac.md#data-processing)
  - [Quality check](atac.md#quality-check)
  - [Trim adapters and filter by base quality](atac.md#trim-adapters-and-filter-by-base-quality)
- [Alignment](atac.md#alignment)
  - [Prepare and index references](atac.md#prepare-and-index-references)



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
- Standard Unix tools
- [cutadapt v1.12](http://cutadapt.readthedocs.io/en/stable/guide.html)




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

Downloaded human hg38 reference genome from [iGenomes](https://emea.support.illumina.com/sequencing/sequencing_software/igenome.html) and HIV genome from the [genomes folder]

hg38 + virus:

```bash
srun --mem 16G --pty /usr/bin/bash

cd /scratchb/sblab/martin03/repository/20200412_liane/reference/fasta

cp /scratcha/sblab/martin03/reference_data/genomes/iGenomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa .
mv genome.fa hg38.fa
cat hg38.fa HIV1_1LTR_circle.fa > hg38_HIV1_1LTR_circle.fa

bwa index hg38_HIV1_1LTR_circle.fa
samtools faidx hg38_HIV1_1LTR_circle.fa

rm hg38.fa

exit
```



### Align and sort

```bash
# Liane_Dupont_SOUK005808
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808/fastq_trimmed_atac

mkdir ../bam_atac

ref=/scratchb/sblab/martin03/repository/20200412_liane/reference/fasta/hg38_HIV1_1LTR_circle.fa

for fq1 in *_R1_001.fastq.gz
do
  bname=${fq1%_L001_R1_001.fastq.gz}
  fq2=${fq1/_R1_/_R2_}
  #echo $bname, $fq1, $fq2
  sbatch -J $bname -o ../bam_atac/$bname.log --mem 32G --wrap "bwa mem -t 20 -M $ref $fq1 $fq2 | \
  samtools view -@ 20 -b - | \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o ../bam_atac/$bname.tmp.bam -"
done

cd ../bam_atac
tail *.log


# Liane_Dupont_SOUK005808_rerun
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/fastq_trimmed_atac

mkdir ../bam_atac

ref=/scratchb/sblab/martin03/repository/20200412_liane/reference/fasta/hg38_HIV1_1LTR_circle.fa

for fq1 in *_R1_001.fastq.gz
do
  bname=${fq1%_L001_R1_001.fastq.gz}
  fq2=${fq1/_R1_/_R2_}
  #echo $bname, $fq1, $fq2
  sbatch -J $bname -o ../bam_atac/$bname.log --mem 32G --wrap "bwa mem -t 20 -M $ref $fq1 $fq2 | \
  samtools view -@ 20 -b - | \
  samtools sort -@ 20 -T /scratchb/sblab/martin03/tmp/$bname -o ../bam_atac/$bname.tmp.bam -"
done

cd ../bam_atac
tail *.log
```



### Merge Liane_Dupont_SOUK005808 and Liane_Dupont_SOUK005808_rerun

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

for bam in *.tmp.bam
do
  bname=${bam%.tmp.bam}
  sbatch -J $bname -o $bname.merge.log --mem 16G --wrap "samtools merge -@ 20 $bname.merge.tmp.bam $bam ../../Liane_Dupont_SOUK005808/bam_atac/$bam"
done

tail *.merge.log
```



### Mark duplicates, index and flagstat

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

mkdir ../flagstat_atac

for bam in *.merge.tmp.bam
do
  bname=${bam%.merge.tmp.bam}
  sbatch -J $bname -o $bname.markdup.log --mem 16G --wrap "sambamba markdup -t 20 $bam $bname.bam 2> $bname.markdup.txt && \
  samtools flagstat $bname.bam > ../flagstat_atac/$bname.txt"
done

ls -lh *.markdup.log

ls *.markdup.txt | xargs head

ls ../flagstat_atac/*.txt | xargs head -n13
# ==> ../flagstat_atac/2-17viA_S4.txt <==
# 207246332 + 0 in total (QC-passed reads + QC-failed reads)
# 909090 + 0 secondary
# 0 + 0 supplementary
# 33959622 + 0 duplicates
# 204217509 + 0 mapped (98.54% : N/A)
# 206337242 + 0 paired in sequencing
# 103168621 + 0 read1
# 103168621 + 0 read2
# 198303606 + 0 properly paired (96.11% : N/A)
# 202011690 + 0 with itself and mate mapped
# 1296729 + 0 singletons (0.63% : N/A)
# 1440700 + 0 with mate mapped to a different chr
# 666366 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/2-17viB_S5.txt <==
# 186900929 + 0 in total (QC-passed reads + QC-failed reads)
# 876123 + 0 secondary
# 0 + 0 supplementary
# 14831000 + 0 duplicates
# 185318054 + 0 mapped (99.15% : N/A)
# 186024806 + 0 paired in sequencing
# 93012403 + 0 read1
# 93012403 + 0 read2
# 179929264 + 0 properly paired (96.72% : N/A)
# 183842798 + 0 with itself and mate mapped
# 599133 + 0 singletons (0.32% : N/A)
# 1091692 + 0 with mate mapped to a different chr
# 447916 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/ConVLPvi_S6.txt <==
# 141066930 + 0 in total (QC-passed reads + QC-failed reads)
# 505818 + 0 secondary
# 0 + 0 supplementary
# 11851171 + 0 duplicates
# 139730692 + 0 mapped (99.05% : N/A)
# 140561112 + 0 paired in sequencing
# 70280556 + 0 read1
# 70280556 + 0 read2
# 136407868 + 0 properly paired (97.05% : N/A)
# 138650150 + 0 with itself and mate mapped
# 574724 + 0 singletons (0.41% : N/A)
# 770646 + 0 with mate mapped to a different chr
# 345178 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/HBxvi_S8.txt <==
# 86573353 + 0 in total (QC-passed reads + QC-failed reads)
# 686949 + 0 secondary
# 0 + 0 supplementary
# 12361218 + 0 duplicates
# 85158110 + 0 mapped (98.37% : N/A)
# 85886404 + 0 paired in sequencing
# 42943202 + 0 read1
# 42943202 + 0 read2
# 82583420 + 0 properly paired (96.15% : N/A)
# 84122012 + 0 with itself and mate mapped
# 349149 + 0 singletons (0.41% : N/A)
# 708168 + 0 with mate mapped to a different chr
# 402979 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/VprVLPvi_S7.txt <==
# 160112840 + 0 in total (QC-passed reads + QC-failed reads)
# 787580 + 0 secondary
# 0 + 0 supplementary
# 18096139 + 0 duplicates
# 157539241 + 0 mapped (98.39% : N/A)
# 159325260 + 0 paired in sequencing
# 79662630 + 0 read1
# 79662630 + 0 read2
# 153164284 + 0 properly paired (96.13% : N/A)
# 155701948 + 0 with itself and mate mapped
# 1049713 + 0 singletons (0.66% : N/A)
# 1160872 + 0 with mate mapped to a different chr
# 528255 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/WTnovi_S1.txt <==
# 156374826 + 0 in total (QC-passed reads + QC-failed reads)
# 610714 + 0 secondary
# 0 + 0 supplementary
# 16990686 + 0 duplicates
# 154908624 + 0 mapped (99.06% : N/A)
# 155764112 + 0 paired in sequencing
# 77882056 + 0 read1
# 77882056 + 0 read2
# 151342374 + 0 properly paired (97.16% : N/A)
# 153782432 + 0 with itself and mate mapped
# 515478 + 0 singletons (0.33% : N/A)
# 807566 + 0 with mate mapped to a different chr
# 352627 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/WTviA_S2.txt <==
# 148512833 + 0 in total (QC-passed reads + QC-failed reads)
# 518449 + 0 secondary
# 0 + 0 supplementary
# 18238030 + 0 duplicates
# 147018579 + 0 mapped (98.99% : N/A)
# 147994384 + 0 paired in sequencing
# 73997192 + 0 read1
# 73997192 + 0 read2
# 144133906 + 0 properly paired (97.39% : N/A)
# 145877676 + 0 with itself and mate mapped
# 622454 + 0 singletons (0.42% : N/A)
# 742214 + 0 with mate mapped to a different chr
# 341813 + 0 with mate mapped to a different chr (mapQ>=5)
#
# ==> ../flagstat_atac/WTviB_S3.txt <==
# 139662646 + 0 in total (QC-passed reads + QC-failed reads)
# 530378 + 0 secondary
# 0 + 0 supplementary
# 11511887 + 0 duplicates
# 138380552 + 0 mapped (99.08% : N/A)
# 139132268 + 0 paired in sequencing
# 69566134 + 0 read1
# 69566134 + 0 read2
# 134736678 + 0 properly paired (96.84% : N/A)
# 137344246 + 0 with itself and mate mapped
# 505928 + 0 singletons (0.36% : N/A)
# 765626 + 0 with mate mapped to a different chr
# 339258 + 0 with mate mapped to a different chr (mapQ>=5)
```

id | mapping rate (%) | duplication rate (%)
:---: | :---: | :---:
WTnovi | 99 | 100*(16990686/154908624) = 11
WTviA | 99 | 100*(18238030/147018579) = 12
WTviB | 99 | 100*(11511887/138380552) = 8
2-17viA | 99 | 100*(33959622/204217509) = 17
2-17viB | 99 | 100*(14831000/185318054) = 8
ConVLPvi | 99 | 100*(11851171/139730692) = 8
VprVLPvi | 98 | 100*(18096139/157539241) = 11
HBxvi | 98 | 100*(12361218/85158110) = 15



### Filter and index

First, create whitelist from blacklist:

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/reference/annotation
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
pigz -d hg38-blacklist.v2.bed.gz
cut -f1-2 ../fasta/hg38_HIV1_1LTR_circle.fa.fai | awk -v OFS="\t" '{print $1, 0, $2}' | bedtools subtract -a - -b hg38-blacklist.v2.bed | sort -k1,1 -k2,2n > hg38-whitelist.v2.bed
grep -v "chrM" hg38-whitelist.v2.bed > hg38-whitelist.v2.nochrM.bed
```

Removing duplicate and supplementary alignments, reads mapping to chrM and blacklisted regions as well as filtering by mapping quality too:

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

for bam in `ls *.bam | grep -v "tmp"`
do
  bname=${bam%.bam}
  sbatch -J $bname -o $bname.clean.log --mem 16G --wrap "samtools view -b -L /scratchb/sblab/martin03/repository/20200412_liane/reference/annotation/hg38-whitelist.v2.nochrM.bed -q 10 -F 3844 -@ 20 $bam > $bname.clean.bam && \
  samtools index $bname.clean.bam"
done

ls -lh *.clean.log

for bam in `ls *.bam | grep -v "tmp"`
do
  n=`samtools view -@ 20 $bam -c`
  echo $bam, $n
done
# 2-17viA_S4.bam, 207246332
# 2-17viA_S4.clean.bam, 150203529
# 2-17viB_S5.bam, 186900929
# 2-17viB_S5.clean.bam, 147706403
# ConVLPvi_S6.bam, 141066930
# ConVLPvi_S6.clean.bam, 113753238
# HBxvi_S8.bam, 86573353
# HBxvi_S8.clean.bam, 62599174
# VprVLPvi_S7.bam, 160112840
# VprVLPvi_S7.clean.bam, 122200958
# WTnovi_S1.bam, 156374826
# WTnovi_S1.clean.bam, 122393380
# WTviA_S2.bam, 148512833
# WTviA_S2.clean.bam, 114883391
# WTviB_S3.bam, 139662646
# WTviB_S3.clean.bam, 112144618
```



### chrM abundance

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

for bam in `ls *.bam | grep -v "clean"`
do
  n_tot=`samtools view -@ 20 $bam -c`
  n_chrM=`samtools view -@ 20 $bam chrM -c`
  pct=`echo "scale=1; 100*$n_chrM/$n_tot" | bc`
  echo -e "$bam\t$n_tot\t$n_chrM\t$pct%"
done
# 2-17viA_S4.bam	207246332	8748662	4.2%
# 2-17viB_S5.bam	186900929	8482179	4.5%
# ConVLPvi_S6.bam	141066930	6728839	4.7%
# HBxvi_S8.bam	86573353	7275949	8.4%
# VprVLPvi_S7.bam	160112840	9867980	6.1%
# WTnovi_S1.bam	156374826	7596478	4.8%
# WTviA_S2.bam	148512833	7975936	5.3%
# WTviB_S3.bam	139662646	7431727	5.3%
```



### Number of human and virus reads, and depth

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

for bam in *.clean.bam
do
  n_human=`samtools idxstats $bam | grep chr* | cut -f1 | xargs samtools view -@ 20 $bam -c`
  depth_human=`samtools depth -r chr1 $bam | awk '{sum+=$3} END {print sum/NR}'`
  n_virus=`samtools view -@ 20 $bam HIV1_1LTR_circle -c`
  depth_virus=`samtools depth -r HIV1_1LTR_circle $bam | awk '{sum+=$3} END {print sum/NR}'`
  echo -e "$bam\t$n_human\t$depth_human\t$n_virus\t$depth_virus"
done
# 2-17viA_S4.clean.bam	150149815	6.09208	53714	1174.05
# 2-17viB_S5.clean.bam	147679650	5.62073	26753	550.877
# ConVLPvi_S6.clean.bam	113731863	5.00402	21375	443.708
# HBxvi_S8.clean.bam	62494430	3.61641	104744	2685.71
# VprVLPvi_S7.clean.bam	122157863	5.43452	43095	973.657
# WTnovi_S1.clean.bam	122393024	5.10915	356	6.48898
# WTviA_S2.clean.bam	114852578	5.12612	30813	708.548
# WTviB_S3.clean.bam	112124981	4.82495	19637	408.762
```



### Insert size plots

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

mkdir ../stats_atac

for bam in *.clean.bam
do
  bname=${bam%.bam}
  mkdir ../stats_atac/$bname
  sbatch -J $bname -o ../stats_atac/$bname/$bname.log --mem 8G --wrap "samtools stats $bam > ../stats_atac/$bname/$bname.txt && \
  ~/sw/samtools/samtools-1.3.1/bin/plot-bamstats -p ../stats_atac/$bname/ ../stats_atac/$bname/$bname.txt"
done

cd ../stats_atac
find . -name "*.clean.txt" | xargs grep "insert size average"
# ./VprVLPvi_S7.clean/VprVLPvi_S7.clean.txt:SN	insert size average:	232.6
# ./ConVLPvi_S6.clean/ConVLPvi_S6.clean.txt:SN	insert size average:	243.5
# ./2-17viB_S5.clean/2-17viB_S5.clean.txt:SN	insert size average:	243.2
# ./2-17viA_S4.clean/2-17viA_S4.clean.txt:SN	insert size average:	231.0
# ./WTviA_S2.clean/WTviA_S2.clean.txt:SN	insert size average:	224.6
# ./HBxvi_S8.clean/HBxvi_S8.clean.txt:SN	insert size average:	248.7
# ./WTviB_S3.clean/WTviB_S3.clean.txt:SN	insert size average:	253.8
# ./WTnovi_S1.clean/WTnovi_S1.clean.txt:SN	insert size average:	231.7
```



### Create bigwig files

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

mkdir ../bw_atac

for bam in *.clean.bam
do
  bname=${bam%.bam}
  sbatch -J $bname -o ../bw_atac/$bname.log --mem 32G --wrap "bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --binSize 1 -p 20 --normalizeUsing CPM"
done

tail ../bw_atac/*.log
```



### Create bigwig files with qPCR fold-change normalisation

Extracting HIV1_1LTR_circle from 2-17viA, WTviA, VprVLPvi and ConVLPvi first:

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

for bam in {2-17viA_S4,WTviA_S2,VprVLPvi_S7,ConVLPvi_S6}.clean.bam
do
  bname=${bam%.bam}
  echo $bam
  samtools view -b -@ 20 $bam HIV1_1LTR_circle > ${bname}.HIV1_1LTR_circle.bam && samtools index ${bname}.HIV1_1LTR_circle.bam  
done
```

Obtain bigwig files using the `--scaleFactor` option:

```bash
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac

# 2-17viA
bam=2-17viA_S4.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1.17 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# WTviA
bam=WTviA_S2.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# VprVLPvi
bam=VprVLPvi_S7.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &

# ConVLPvi
bam=ConVLPvi_S6.clean.HIV1_1LTR_circle.bam
bname=${bam%.bam}
nohup bamCoverage -b $bam -o ../bw_atac/$bname.bw -of bigwig --scaleFactor 1.5 --binSize 1 -p 20 --normalizeUsing CPM > ../bw_atac/$bname.log &
```



### Clean-up

```bash
# Liane_Dupont_SOUK005808
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808
rm -r fastq_trimmed_atac bam_atac

# Liane_Dupont_SOUK005808_rerun
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun
rm -r fastq_trimmed_atac
cd bam_atac
rm *.tmp.bam
```



### ATACseqQC

Interested in the fragment size distribution plots:

```r
#srun --mem 16G --pty /usr/bin/bash
#cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/bam_atac
#/home/martin03/sw/R/R-3.5.2/bin/R

# install libs
#library(BiocManager)
#BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "BSgenome.Hsapiens.UCSC.hg38", "TxDb.Hsapiens.UCSC.hg38.knownGene", "phastCons100way.UCSC.hg38"))

# load lib
library(ATACseqQC)

# estimate library complexity
f <- "WTnovi_S1.clean.bam"
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
