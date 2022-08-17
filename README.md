# Vamos: VNTR Annotation tool using efficient MOtifs Sets

## Updates 

## Getting Started

To install vamos, g++ (>= 8.3.0), htslib, abpoa, edlib and alglib are required.
htslib and abpoa can be installed through bioconda.
Static libraries libalglib.a and libedlib.a are distributed along with vamos.

Install required libraries through conda
```
conda create --name vamos python=3.10
conda activate vamos
conda install -c bioconda --file requirements.txt
```

Download the [latest release](https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.0.0.tar.gz)
```
wget https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.0.0.tar.gz
tar -zxvf vamos-v1.0.0.tar.gz && cd vamos-v1.0.0
```

Or download the latest code from github
```
git clone https://github.com/ChaissonLab/vamos.git
```

Make from source and run with test data:
```
cd vamos*/src/ && make
vamos --readwise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/readwise.bed -t 16
vamos --locuswise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/locuswise.bed -t 16
vamos --single_seq -b ../example/one_read.fasta -r ../example/one_region_motif.bed -s NA24385_CCS_h1 -o ../example/single_seq.vcf (ONLY SUPPORT SINGLE THREAD!)
```

## Table of Contents

- [Introduction](#introduction)
  - [Efficient motif set](#emotif)
- [Installation](#install)
  - [Building vamos from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate the annotation for each read](#readwise)
  - [To generate the consensus annotation for each VNTR locus](#locuswise)
  - [To generate annotation for one single read](#single_seq)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [BED for readwise mode](#BED)
  - [VCF for locuswise mode](#VCF)

## <a name="introduction"></a>Introduction
Vamos is a tool to perform run-length encoding of VNTR sequences using a set of selected motifs from all motifs observed at that locus.
Vamos guarantees that the encoding sequence is winthin a bounded edit distance of the original sequence. 
For example, a VNTR sequence ACGGTACTGTACGT may be encoded to a more compact representation - ACGGT|AC**G**GT|ACGT using motif set [ACGGT, ACGT].

Vamos can generate annotation for each read, given a set of motifs for each VNTR locus. (`--readwise` mode)
Vamos can generate annotation for each VNTR locus by aggragating annotations of reads. (`--locuswise` mode)
Vamos can generate annotation for subsequence of a read lifted from a particular VNTR locus. (`--single_seq` mode)

### <a name="emotif"></a>Efficient motif set
We generate a set of efficient VNTR motifs from 32 haplotype-resolved LRS genomes sequenced for population references and diversity panels.
XXXX (some statistics about # of VNTRs)
XXXX (path)

## <a name="install"></a>Installation
### <a name="build"></a>Building vamos from source files
You can build vamos from source files. 
Make sure to install the required libraries: g++ (>= 8.3.0), htslib, abpoa, edlib and alglib before compiling. 

Install required libraries
```
conda create --name vamos python=3.10
conda activate vamos
conda install -c bioconda --file requirements.txt
```

Download the [latest release](https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.0.0.tar.gz)
```
wget https://github.com/ChaissonLab/vamos/archive/refs/tags/vamos-v1.0.0.tar.gz
tar -zxvf vamos-v1.0.0.tar.gz
cd vamos-v1.0.0/src; make
```
Or, you can use `git clone` command to download the source code.
This gives you the latest version of vamos, which might be still under development.
```
git clone https://github.com/ChaissonLab/vamos.git
cd vamos/src; make
```

### <a name="binary"></a>Pre-built binary executable file for Linux/Unix 
If you meet any compiling issue, please try the pre-built binary file:
XXXX(TODO)
```
wget https://github.com/ChaissonLab/vamos/releases/download/v1.0.0/vamos-v1.0.0_x64-linux.tar.gz
tar -zxvf vamos-v1.0.0_x64-linux.tar.gz
```

## <a name="usage"></a>General usage
### <a name="readwise"></a>To generate the annotation for each read
```
vamos --readwise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/readwise.bed -t 16
```
### <a name="locuswise"></a>To generate the consensus annotation for each VNTR locus
```
vamos --locuswise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/locuswise.bed -t 16
```
### <a name="single_seq"></a>To generate annotation for one single read
```
vamos --single_seq -b ../example/one_read.fasta -r ../example/one_region_motif.bed -s NA24385_CCS_h1 -o ../example/single_seq.vcf 
```

## <a name="cmd"></a>Commands and options
```
Usage: vamos [subcommand] [options] [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf/bed] [-s sample_name] [-t threads] 
Version: v1.0.0
subcommand:
vamos --readwise   [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.bed] [-s sample_name] [-t threads] 
vamos --locuswise  [-b in.bam] [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] [-t threads] 
vamos --single_seq [-b in.fa]  [-r vntrs_region_motifs.bed] [-o output.vcf] [-s sample_name] (ONLY FOR SINGLE LOCUS!!) 
   Input: 
       -b   FILE         input indexed bam file. 
       -r   FILE         file containing region coordinate and motifs of each VNTR locus. 
                         The file format: columns `chrom,start,end,motifs` are tab-delimited. 
                         Column `motifs` is a comma-separated (no spaces) list of motifs for this VNTR. 
       -s   CHAR         sample name. 
   Output: 
       -o   FILE         output bed/vcf file. 
   Dynamic Programming: 
       -d   DOUBLE       penalty of indel in dynamic programming (double) DEFAULT: 1.0. 
       -c   DOUBLE       penalty of mismatch in dynamic programming (double) DEFAULT: 1.0. 
       -a   DOUBLE       Global accuracy of the reads. DEFAULT: 0.98. 
       --naive           specify the naive version of code to do the annotation, DEFAULT: faster implementation. 
   Aggregate Annotation: 
       -f   DOUBLE       filter noisy read annotations, DEFAULT: 0.0 (no filter). 
       --clust           use hierarchical clustering to judge if a VNTR locus is het or hom. 
   Others: 
       -t   INT          number of threads, DEFAULT: 1. 
       --debug           print out debug information. 
       -h                print out help message. 
```

## <a name="input"></a>Input

## <a name="output"></a>Output
### <a name="BED"></a>BED for readwise mode
### <a name="VCF"></a>VCF for locuswise mode
