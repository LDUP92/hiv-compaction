
## Contents

- [Libraries][atac.md#libraries]
- [Software](atac.md#software-and-operating-system)
- [Processing](atac.md#processing)
  - [Quality check](atac.md#quality-check)



## Libraries

id | files
---|------
WTnoVi | WTnoVi.R1.fastq.gz,WTnoVi.R2.fastq.gz
WTnoVLP | WTnoVLP.R1.fastq.gz,WTnoVLP.R2.fastq.gz
WTConVLP | WTConVLP.R1.fastq.gz,WTConVLP.R2.fastq.gz
WTVprVLP | WTVprVLP.R1.fastq.gz,WTVprVLP.R2.fastq.gz



## Software and operating system

Software:

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)



## Processing

### Quality check

```bash
# Liane_Dupont_SOUK005808
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808/fastq_atac

mkdir ../fastqc_atac

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc_atac/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc_atac $fq"
done


# Liane_Dupont_SOUK005808_rerun
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/fastq_atac

mkdir ../fastqc_atac

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc_atac/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc_atac $fq"
done
```

In `C02Q70MUFVH8`,

```bash
cd /Users/martin03/Desktop/

# Liane_Dupont_SOUK005808
mkdir Liane_Dupont_SOUK005808/ && cd Liane_Dupont_SOUK005808/
rsync -arvuP martin03@10.20.236.34:/scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808/fastqc_atac/*.html .

# Liane_Dupont_SOUK005808_rerun
cd ..
mkdir Liane_Dupont_SOUK005808_rerun/ && cd Liane_Dupont_SOUK005808_rerun/
rsync -arvuP martin03@10.20.236.34:/scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/fastqc_atac/*.html .
```

- The quality looks fine overall
- Nextera adapters
- High concentration of Gs > 100bp
- Interesting pattern at < 15bp



### Trim adapters and filter by base quality

```bash
# Liane_Dupont_SOUK005808
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808/fastq_atac

mkdir ../fastq_trimmed_atac

for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed_atac/$bname.log --mem 4G --wrap "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 -o ../fastq_trimmed_atac/$fq1 -p ../fastq_trimmed_atac/$fq2 $fq1 $fq2 > ../fastq_trimmed_atac/$bname.txt"
done

cd ../fastq_trimmed_atac

ls -lh *.log

grep "Read 1 with adapter" *.txt
# 2-17viA_S4.txt:  Read 1 with adapter:              26,910,433 (66.4%)
# 2-17viB_S5.txt:  Read 1 with adapter:              26,019,657 (76.5%)
# ConVLPvi_S6.txt:  Read 1 with adapter:              15,134,338 (65.0%)
# HBxvi_S8.txt:  Read 1 with adapter:               8,264,339 (59.7%)
# VprVLPvi_S7.txt:  Read 1 with adapter:              20,441,839 (67.4%)
# WTnovi_S1.txt:  Read 1 with adapter:              15,727,983 (67.2%)
# WTviA_S2.txt:  Read 1 with adapter:              14,829,715 (61.3%)
# WTviB_S3.txt:  Read 1 with adapter:              14,710,570 (66.2%)

grep "Read 2 with adapter" *.txt
# 2-17viA_S4.txt:  Read 2 with adapter:              25,796,060 (63.7%)
# 2-17viB_S5.txt:  Read 2 with adapter:              26,391,888 (77.6%)
# ConVLPvi_S6.txt:  Read 2 with adapter:              14,692,925 (63.1%)
# HBxvi_S8.txt:  Read 2 with adapter:               8,129,417 (58.8%)
# VprVLPvi_S7.txt:  Read 2 with adapter:              19,379,084 (63.9%)
# WTnovi_S1.txt:  Read 2 with adapter:              16,218,370 (69.3%)
# WTviA_S2.txt:  Read 2 with adapter:              14,354,102 (59.3%)
# WTviB_S3.txt:  Read 2 with adapter:              14,703,540 (66.1%)

grep "Pairs written (passing filters):" *.txt
# 2-17viA_S4.txt:Pairs written (passing filters):    40,495,353 (100.0%)
# 2-17viB_S5.txt:Pairs written (passing filters):    34,011,435 (100.0%)
# ConVLPvi_S6.txt:Pairs written (passing filters):    23,278,368 (100.0%)
# HBxvi_S8.txt:Pairs written (passing filters):    13,830,003 (100.0%)
# VprVLPvi_S7.txt:Pairs written (passing filters):    30,317,993 (100.0%)
# WTnovi_S1.txt:Pairs written (passing filters):    23,411,106 (100.0%)
# WTviA_S2.txt:Pairs written (passing filters):    24,188,470 (100.0%)
# WTviB_S3.txt:Pairs written (passing filters):    22,229,914 (100.0%)



# Liane_Dupont_SOUK005808_rerun
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun/fastq_atac

mkdir ../fastq_trimmed_atac

for fq1 in *R1*.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%_L001_R1_001.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed_atac/$bname.log --mem 4G --wrap "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 -o ../fastq_trimmed_atac/$fq1 -p ../fastq_trimmed_atac/$fq2 $fq1 $fq2 > ../fastq_trimmed_atac/$bname.txt"
done

cd ../fastq_trimmed_atac

ls -lh *.log

grep "Read 1 with adapter" *.txt
# 2-17viA_S4.txt:  Read 1 with adapter:              34,574,212 (55.2%)
# 2-17viB_S5.txt:  Read 1 with adapter:              40,673,727 (68.9%)
# ConVLPvi_S6.txt:  Read 1 with adapter:              24,723,522 (52.6%)
# HBxvi_S8.txt:  Read 1 with adapter:              13,800,507 (47.4%)
# VprVLPvi_S7.txt:  Read 1 with adapter:              27,706,945 (56.1%)
# WTnovi_S1.txt:  Read 1 with adapter:              32,536,129 (59.7%)
# WTviA_S2.txt:  Read 1 with adapter:              25,259,002 (50.7%)
# WTviB_S3.txt:  Read 1 with adapter:              25,888,539 (54.7%)

grep "Read 2 with adapter" *.txt
# 2-17viA_S4.txt:  Read 2 with adapter:              34,350,538 (54.8%)
# 2-17viB_S5.txt:  Read 2 with adapter:              41,051,831 (69.6%)
# ConVLPvi_S6.txt:  Read 2 with adapter:              24,656,843 (52.5%)
# HBxvi_S8.txt:  Read 2 with adapter:              13,831,921 (47.5%)
# VprVLPvi_S7.txt:  Read 2 with adapter:              27,435,041 (55.6%)
# WTnovi_S1.txt:  Read 2 with adapter:              33,015,436 (60.6%)
# WTviA_S2.txt:  Read 2 with adapter:              25,184,872 (50.6%)
# WTviB_S3.txt:  Read 2 with adapter:              26,050,430 (55.0%)

grep "Pairs written (passing filters):" *.txt
# 2-17viA_S4.txt:Pairs written (passing filters):    62,673,268 (100.0%)
# 2-17viB_S5.txt:Pairs written (passing filters):    59,000,968 (100.0%)
# ConVLPvi_S6.txt:Pairs written (passing filters):    47,002,188 (100.0%)
# HBxvi_S8.txt:Pairs written (passing filters):    29,113,199 (100.0%)
# VprVLPvi_S7.txt:Pairs written (passing filters):    49,344,637 (100.0%)
# WTnovi_S1.txt:Pairs written (passing filters):    54,470,950 (100.0%)
# WTviA_S2.txt:Pairs written (passing filters):    49,808,722 (100.0%)
# WTviB_S3.txt:Pairs written (passing filters):    47,336,220 (100.0%)
```



### Clean-up

```bash
# Liane_Dupont_SOUK005808
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808
rm -r fastq_atac fastqc_atac

# Liane_Dupont_SOUK005808_rerun
cd /scratchb/sblab/martin03/repository/20200412_liane/data/20200422/Liane_Dupont_SOUK005808_rerun
rm -r fastq_atac fastqc_atac
```
