# Vamos: VNTR Annotation tool using efficient MOtifs Sets

## Updates 

## Getting Started

To install vamos, g++ (>= 8.3.0), htslib, abpoa, edlib and alglib are required.
htslib and abpoa can be installed by bioconda.
Static libraries libalglib.a and libedlib.a are distributed along with vamos.

Download the [latest release] (https://github.com/ChaissonLab/vamos/XXX)
```
wget https://github.com/ChaissonLab/vamos.git/releases/download/XXXX
tar -zxvf xxxx.tar.gz && cd xxxx
```

Or download the latest code from github
```
git clone https://github.com/ChaissonLab/vamos.git
```

Install required libraries
```
conda create --name vamos python=3.10
conda activate vamos
conda install -c bioconda --file requirements.txt
```

Make from source and run with test data:
```
cd vamos/src/ && make
vamos --readwise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/readwise.bed -t 16
vamos --locuswise -b ../example/toy.bam -r ../example/region_motif.bed -s NA24385_CCS_h1 -o ../example/locuswise.bed -t 16
vamos --single_seq -b ../example/one_read.fasta -r ../example/one_region_motif.bed -s NA24385_CCS_h1 -o ../example/single_seq.vcf 
```

## Table of Contents

- [Introduction](#introduction)
- [Installation](#install)
  - [Building vamos from source files](#build)
  - [Pre-built binary executable file for Linux/Unix](#binary)
- [General usage](#usage)
  - [To generate the annotation for each read](#readwise)
  - [To generate the consensus annotation for each VNTR locus](#locuswise)
  - [To generate annotation for one single read] (#single_seq)
- [Commands and options](#cmd)
- [Input](#input)
- [Output](#output)
  - [BED for readwise mode](#BED)
  - [VCF for locuswise mode](#VCF)

## <a name="introduction"></a>Introduction
Vamos is a tool to perform run-length encoding of VNTR sequences using a set of selected motifs from all motifs observed at that locus.
Vamos gauranttes that the encoded sequence (annotation) is winthin a bounded edit distance of the original sequence. 
For example, a VNTR sequence `ACGGTACTGTACGT` may be encoded to a more compact representation `ACGGT, AC\textcolor{red}{G}GT, ACGT` using motif set `[ACGGT, ACGT]`

Vamos can generate annotation for each read at a VNTR locus, given a set of motifs at that locus. (`--readwise` mode)
Vamos can generate annotation for each VNTR locus by aggragating annotations of reads. (`--locuswise` mode)
Vamos can generate annotation for subsequence of a read lifted from a VNTR locus. (`--single_seq` mode)

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

Download the [latest release] (https://github.com/ChaissonLab/vamos/XXX)
```
wget https://github.com/xxxx
tar -zxvf vamos-xxx.tar.gz
cd vamos/src; make

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
XXXX(TODO) put up the help page

## <a name="input"></a>Input

## <a name="output"></a>Output
### <a name="BED"></a>BED for readwise mode
### <a name="VCF"></a>VCF for locuswise mode
