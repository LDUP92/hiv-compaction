
## Contents

- [Libraries](screen.md#libraries)
- [Software](screen.md#software)
- [Clip and trim](screen.md#clip-and-trim) 
- [Alignment](screen.md#alignment)
- [MAGeCK](screen.md#mageck)



## Libraries

id | files
---|------
Vpr_LIB | Vpr_LIB.fastq.gz
Vpr_SORT | Vpr_SORT.fastq.gz



## Software

Software:

- [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- [samtools](http://samtools.sourceforge.net/)
- [MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/)



## Clip and trim

```bash
cd ~

fastx_clipper -Q33 -l 15 -n -a GTTTAAGAGCTAAG -v -i Vpr_LIB.fastq -o Vpr_LIB_clipped.fastq
fastx_trimmer -l 19 -v -i Vpr_LIB_clipped.fastq -o Vpr_LIB_trimmed.fastq -Q 33

fastx_clipper -Q33 -l 15 -n -a GTTTAAGAGCTAAG -v -i Vpr_SORT.fq -o Vpr_SORT_clipped.fastq
fastx_trimmer -l 19 -v -i Vpr_SORT_clipped.fastq -o Vpr_SORT_trimmed.fastq -Q 33
```



## Alignment

```bash
cd ~

bowtie2 -x Vprlib -U Vpr_LIB_trimmed.fastq --norc | samtools view -bS - > Vpr_LIB.bam
bowtie2 -x Vprlib -U Vpr_SORT_trimmed.fastq --norc | samtools view -bS - > Vpr_SORT.bam
```



## MAGeCK

```bash
cd ~

mageck count -l Vprlib.csv -n escneg --sample-label "Lib,sort" --fastq Vpr_LIB.bam Vpr_SORT.bam
mageck test -n Vprlib -k escneg.count.txt -t 1 -c 0
```
