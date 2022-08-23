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
For example, a VNTR sequence ACGGT|ACTGT|ACGT may be encoded to a more compact representation: ACGGT|AC**G**GT|ACGT using efficient motif set [ACGGT, ACGT].
The edit distance between the original VNTR sequence and encoding sequence is 1. 

Vamos can generate annotation for each read, given a set of motifs for each VNTR locus. (`--readwise` mode)
Vamos can generate annotation for each VNTR locus by aggragating annotations of reads. (`--locuswise` mode)
Vamos can generate annotation for subsequence of a read lifted from a particular VNTR locus. (`--single_seq` mode)

### <a name="emotif"></a>Efficient motif set
We defined VNTR loci and motifs using a collection of 32 haplotype-resolved LRS genomes constructed by Human Genome Structural Variation Consortium.
XXXX(TOADD) loci of simple repeating sequences on the GRCh38 assembly were obtained from the table browser tool of the UCSC Genome Browser.
For each assembly, VNTR sequences were lifted-over and decomposed into motifs by Tandem Repeats Finder (TRF). Post-filtering step leaves 467104 well-resolved VNTR loci. 
We propose efficient motif set as a smallest set of motifs, such that the string decompositions of the assembly alleles are bounded by a given edit distance. 
[snakefile/configs/vntr_region_motifs.e.bed.gz](https://github.com/ChaissonLab/vamos/blob/master/snakefile/configs/vntr_region_motifs.e.bed.gz) provides efficent motifs for 467104 VNTR loci.
[snakefile/configs/vntr_region_motifs.o.bed.gz](https://github.com/ChaissonLab/vamos/blob/master/snakefile/configs/vntr_region_motifs.o.bed.gz) provides original motifs for 467104 VNTR loci.

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
```
wget https://github.com/ChaissonLab/vamos/releases/download/vamos-v1.0.0/vamos-v1.0.0_x64-linux.tar.gz
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
Vamos works with indexed bam file under `--readwise` and `--locuswise` modes, and works with fasta file under `--single_seq` mode. 

## <a name="output"></a>Output
Vamos outputs BED file under `--readwise` mode, and outputs vcf file under `--locuswise` and `--single_seq` modes. 


### <a name="BED"></a>BED for readwise mode
Vamos generates annotation for each read in BED file under `--readwise` mode.

The BED file contains five columns - `CHROM`, `START`, `END`, `MOTIFS`, `INFO`.
`CHROM:START-END`:coordinate of VNTR locus. \\
`MOTIFS`: a comma-separated list of motifs at the locus. \\ 
`INFO`: a colon-separated list of information about the read annotation. (Read_name:Haplotype(0 - not determined, 1/2 - haplotype):Annotation_length:Annotation(comma-separated, no spaces):Read_lifted_seq) \\

The following is an example of the BED file: 
```
##fileformat=BED
##source=vamos_v1.0.0
##INFO=<Read_name:Haplotype(0 - not determined, 1/2 - haplotype):Annotation_length:Annotation(comma-separated, no spaces):Read_lifted_seq;,Description="read annotation information per read">
#CHROM	START	END	MOTIFS	INFO
chr1	189828	189966	TGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAGGGCAGGAGGAGGGTGTGGGATGGTGGAGGGGTT,TGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAAGGGAGGGGGAGGATGTGGGATGGTGGAGGGG,GAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAAGGGAGGGGGAGGATGTGGGATGGTGGAGGGGGA	cluster2_000015F:0:2:0,1:ATGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAGGGCAGGAGGAGGGTGTGGGATGGTGGAGGGGTTTGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAAGGGAGGGGGAGGATGTGGGATGGTGGAGGGG;cluster2_000015F:0:2:0,1:ATGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAGGGCAGGAGGAGGGTGTGGGATGGTGGAGGGGTTTGAGAAGGCAGAGGCGCGACTGGGGTTCATGAGGAAAGGGAGGGGGAGGATGTGGGATGGTGGAGGGG;
```

### <a name="VCF"></a>VCF for locuswise and single_seq modes
Vamos generates annotation for each VNTR locus in VCF file under `--locuswise` and `--single_seq` modes.


The following shows an example of the VCF file
```
##fileformat=VCFv4.1
##source=vamos_V1.0.0
##contig=<ID=chr1,length=248956422>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">
##INFO=<ID=RU,Number=1,Type=String,Description="Comma separated motif sequences list in the reference orientation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=ALTANNO_H1,Number=1,Type=String,Description=""Motif representation for the h1 alternate allele>
##INFO=<ID=ALTANNO_H2,Number=1,Type=String,Description=""Motif representation for the h2 alternate allele>
##INFO=<ID=LEN_H1,Number=1,Type=Integer,Description=""Length of the motif annotation for the h1 alternate allele>
##INFO=<ID=LEN_H2,Number=1,Type=Integer,Description=""Length of the motif annotation for the h2 alternate allele>
##FILTER=<ID=PASS,Description="All filters passed">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##ALT=<ID=VNTR,Description="Allele comprised of VNTR repeat units">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA24385_CCS_h1
chr1	191351	.	N	<VNTR>	.	PASS	END=191386;RU=CACCACAGAAAACAGAG,CACCACAGAAAACAGAGC;SVTYPE=VNTR;ALTANNO_H1=1,0;LEN_H1=2;	GT	1/1
```
