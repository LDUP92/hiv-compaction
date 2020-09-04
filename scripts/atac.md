
## Contents

- [Libraries](atac.md#libraries)
- [Software and operating system](atac.md#software-and-operating-system)
- [Data processing](atac.md#data-processing)
  - [Quality check](atac.md#quality-check)
  - [Trim adapters and filter by base quality](atac.md#trim-adapters-and-filter-by-base-quality)



## Libraries

id | files
---|------
WTnoVi | WTnoVi.R1.fastq.gz,WTnoVi.R2.fastq.gz
WTnoVLP | WTnoVLP.R1.fastq.gz,WTnoVLP.R2.fastq.gz
WTConVLP | WTConVLP.R1.fastq.gz,WTConVLP.R2.fastq.gz
WTVprVLP | WTVprVLP.R1.fastq.gz,WTVprVLP.R2.fastq.gz

Stored in folder `~/fastq_atac`



## Software and operating system

Software:

- [FastQC v0.11.3](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Standard Unix tools



Operating system:

  - CentOS Linux release 7.3.1611 (OS used during code development)
  - [slurm job scheduling system v19.05.0](https://slurm.schedmd.com/quickstart.html)



## Data processing

### Quality check

```bash
cd fastq_atac

mkdir ../fastqc_atac

for fq in *.fastq.gz
do
  bname=${fq%.fastq.gz}
  sbatch -J $bname -o ../fastqc_atac/$bname.log --mem 4G --wrap "fastqc --noextract --nogroup -q -o ../fastqc_atac $fq"
done
```


### Trim adapters and filter by base quality

```bash
cd fastq_atac

mkdir ../fastq_trimmed_atac

for fq1 in *.R1.fastq.gz
do
  fq2=${fq1/_R1_/_R2_}
  bname=${fq1%.R1.fastq.gz}
  sbatch -J $bname -o ../fastq_trimmed_atac/$bname.log --mem 4G --wrap "cutadapt -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -m 15 -q 20 -o ../fastq_trimmed_atac/$fq1 -p ../fastq_trimmed_atac/$fq2 $fq1 $fq2 > ../fastq_trimmed_atac/$bname.txt"
done
```
